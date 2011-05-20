/*
TinyTinyJPG library
Copyright (c) 2011, Rasmus Neckelmann (neckelmann@gmail.com)
All rights reserved.

This source code is released under the terms of the BSD license; see the 
file LICENSE.TXT for details.
*/
#include <math.h>
#if defined(EMMINTRIN)
  #include <emmintrin.h>
#endif

#include "TinyTinyJPG.h"

/* Clamp x to 0-255 */
#define BYTECLAMP(x) ((static_cast<unsigned int>(x) > 255)?(((~x) >> 31) & 0xff):x)

namespace ttjpg {

  /*===========================================================================
  Global stuff
  ===========================================================================*/
  /* MCU layout description */
  struct MCULayoutDesc {
    int nXSampleFactor[3],nYSampleFactor[3]; /* This is used to lookup a layout */
    
    int nDataUnitsPerMCU;     /* Number of DCT blocks in a MCU */
    int nDUComponents[6];     /* Component-index of each data unit */
    int nMCUXSize,nMCUYSize;  /* Size of MCU in pixels */
  };  
  
  /* These MCU layouts hold information needed to map block to 
     the final image. Actual positioning in hardcoded for each layout in the 
     _UpdateMCU() function */
  MCULayoutDesc g_MCULayouts[JPG_NUM_MCU_LAYOUTS] = {    
    1,1,1,1,1,1,3,0,1,2,0,0,0,8,8,  /* JPG_MCU_LAYOUT_YCBCR_444 */    
    2,1,1,1,1,1,4,0,0,1,2,0,0,16,8, /* JPG_MCU_LAYOUT_YCBCR_422 */    
    2,1,1,2,1,1,6,0,0,0,0,1,2,16,16 /* JPG_MCU_LAYOUT_YCBCR_420 */
  };
  
  /* Function for looking up MCU layouts in the above table */        
  int _LookupMCULayout(ComponentInfo *p,int nNumComponents) {
    if(nNumComponents != 3) return -1; /* Must be 3-component */
    
    for(int i=0;i<JPG_NUM_MCU_LAYOUTS;i++) {
      bool bMatch = true;
      for(int j=0;j<3;j++) {
        if(g_MCULayouts[i].nXSampleFactor[j] != p[j].nXSampleFactor ||
           g_MCULayouts[i].nYSampleFactor[j] != p[j].nYSampleFactor) {
          bMatch = false;
          break;
        }
      }
      
      /* This' the one? */
      if(bMatch)
        return i;
    }
    
    throw DecodeError(ERR_UNSUPPORTED_MCU_LAYOUT);
  }

  /*===========================================================================
  Shared input stream functionality
  ===========================================================================*/
  int InputStream::readWord(void) {
    /* Big endian 16-bit word */
    return (readByte() << 8) | readByte();
  }
  
  int InputStream::readMarker(void) {
    /* Read until we get a marker */
    while(!isEndOfStream()) {
      if(readByte() == 0xff) {
        int nMarker = 0;
        do {
          nMarker = readByte();
          if(isEndOfStream()) return 0;
        } while(nMarker == 0xff);
        
        /* If it's an APP0 marker (EXIF, etc) skip it */
        if(nMarker == 0xE1) {
          int nSkipLen = readWord() - 2;
          skip(nSkipLen);
          return readMarker();
        }
        return nMarker;
      }
    }
    
    return 0;
  }
  
  uint32 InputStream::readBits(int nNum) {
    /* Read the next nNum bits, if we hit a 0xFFxx marker take a 
       note of it (embedded marker) */
    if(nNum == 0) return 0;
        
    if(m_nBits < nNum) {
      /* We don't have enough bits buffered, need to read some more */
      int c1 = readByte();
      if(c1 == 0xff) readByte();            
      int c2 = readByte();
      if(c2 == 0xff) {
        m_nEmbeddedMarker = readByte();             
      }

      m_nBitBuf |= c2 << (16 - m_nBits);      
      m_nBitBuf |= c1 << (24 - m_nBits);      
      m_nBits += 16;    
    }

    /* Grab the bits we need from the buffer */
    uint32 n = m_nBitBuf >> (32 - nNum);
    m_nBitBuf <<= nNum;
    m_nBits -= nNum;
    
    return n;
  }
  
  void InputStream::flushBitStream(void) {
    /* Discard any remaining buffered bits */
    m_nBits = 0;
    m_nBitBuf = 0;    
    m_nEmbeddedMarker = 0;
  }  

  /*===========================================================================
  Huffman decoding tree implementation
  ===========================================================================*/  
  void HuffmanTable::buildTree(void) {
    /* Build huffman decoding tree -- create initial tree */
    Tree.grow(3);
    Tree[0].pChildren[0] = &Tree[1];
    Tree[0].pChildren[1] = &Tree[2];
    
    nNextFreeNode = 1;
    nNumFreeNodes = 2;
    
    /* For each level of the tree... */
    for(int i=0;i<nMaxBits;i++) {
      /* Add all symbols at this level... */
      int *pnSymbols = &nTable[i * 256];
      
      for(int j=0;j<nL[i];j++) {
        if(nNumFreeNodes > 0) {
          /* Insert symbol here */
          Tree[nNextFreeNode].nSymbol = pnSymbols[j];
          nNextFreeNode++;
          nNumFreeNodes--;
        }
        else {
          /* Invalid tree... can't possible build it like this */
          bDefined = false;
          throw DecodeError(ERR_INVALID_HUFFMAN_TREE);
        }        
      }
      
      /* Every remaining node at this level will be a splitter */
      int nNextLevelNode = nNextFreeNode + nNumFreeNodes;
      int nNewNodes = 0;
      
      if(nNumFreeNodes > 0) {
        /* Every new splitter will have 2 children... */
        nNewNodes = nNumFreeNodes * 2;
        Tree.grow(nNewNodes);
      }

      while(nNumFreeNodes > 0) {
        Tree[nNextFreeNode].pChildren[0] = &Tree[nNextLevelNode++];
        Tree[nNextFreeNode].pChildren[1] = &Tree[nNextLevelNode++];        
        
        nNextFreeNode++;
        nNumFreeNodes--;        
      }
      
      nNextFreeNode = Tree.size() - nNewNodes;
      nNumFreeNodes = nNewNodes;
    }        
  }
  
  int HuffmanTable::decode(InputStream *pIn,int *pnLeadingZeros) {
    /* Decode the next value. If pnLeadingZeros is NULL, we'll 
       decode a DC value, AC otherwise */
    int nCurrentNode = 0;
    HuffmanNode *pCurrentNode = &Tree[0];
    
    if(pnLeadingZeros != NULL) *pnLeadingZeros = 0;
    
    while(pCurrentNode != NULL) {
      if(pCurrentNode->nSymbol >= 0) {
        /* Leaf node */
        int nL = pCurrentNode->nSymbol;
        
        if(nL == 0) return JPG_HUFFMAN_EOB; /* EOB */
        if(nL == 0xf0) return JPG_HUFFMAN_ZRL; /* 16 zeros */
        
        int nAmp = nL;
        
        if(pnLeadingZeros != NULL) { 
          *pnLeadingZeros = (nAmp & 0xf0) >> 4;
          nAmp &= 0x0f;                
        }
                
        uint32 n = pIn->readBits(nAmp);
        if(nAmp == 1) return n ? 1 : -1;        
        int nS = 1 << (nAmp - 1);        
        if(n & nS) return nS + (n & (nS-1));
        return -(1 << nAmp) + (n & (nS-1)) + 1;
      }
      else {
        /* Read bit to determine which way to go in tree */
        pCurrentNode = pCurrentNode->pChildren[pIn->readBits(1)];
      }
    }
    
    throw DecodeError(ERR_INVALID_HUFFMAN_CODE);
  }

  /*===========================================================================
  JPEG decoder class implementation
  ===========================================================================*/  
  Decoder::Decoder() {
    /* Initialize ZZ-ordering */
    _InitZZOrdering();
    
    /* Make YCBCR->RGB tables */
    _InitRGBTables();
  }
  
  Decoder::~Decoder() {
  }
      
  void Decoder::readImageInfo(InputStream *pIn,RGBImage *pRGB) {
    Info Info;
  
    /* Parse stream just for image info */
    _Parse(pIn,pRGB,&Info,true);
  }
  
  void Decoder::readImage(InputStream *pIn,RGBImage *pRGB) {
    Info Info;

    /* Parse JPG file incl. full image */
    _Parse(pIn,pRGB,&Info,false);
  }
  
  void Decoder::_Parse(InputStream *pIn,RGBImage *pRGB,Info *pi,bool bGetInfoOnly) {
    /* SOI marker must be first */
    if(pIn->readMarker() != 0xD8)
      throw DecodeError(ERR_NOT_JPEG);
      
    /* Read markers until start of frame */
    bool bSOF = false;
    
    while(!pIn->isEndOfStream()) {
      int nMarker = pIn->readMarker();
      
      if(nMarker == 0xC0) {
        /* Start Of Frame */
        _ParseSOF(pIn,pRGB,pi);
        
        pi->nMCUIndex = 0;
        bSOF = true;
        break;
      }
      else if(!bGetInfoOnly) {
        _InterpretMarker(pIn,pRGB,pi,nMarker);
      }
    }
    
    if(!bSOF)
      throw DecodeError(ERR_NO_FRAME);
    
    /* If we're only getting image info stop here */
    if(!bGetInfoOnly) {
      /* Allocate room for image */
      pRGB->release();
      pRGB->pc = new uint8[pRGB->nWidth * pRGB->nHeight * 3];
      memset(pRGB->pc,0,pRGB->nWidth * pRGB->nHeight * 3);
        
      /* Decode frame */
      _DecodeFrame(pIn,pRGB,pi);
    }
  }
  
  void Decoder::_InterpretMarker(InputStream *pIn,RGBImage *pRGB,Info *pi,int nMarker) {
    switch(nMarker) {
      case 0xC4: /* Define Huffman Tables */
        _ParseDHT(pIn,pRGB,pi);
        break;
      case 0xDB: /* Define Quantization Tables */
        _ParseDQT(pIn,pRGB,pi);
        break;
      case 0xDD: /* Define Restart Interval */
        _ParseDRI(pIn,pRGB,pi);
        break;
    }
  }
  
  void Decoder::_DecodeFrame(InputStream *pIn,RGBImage *pRGB,Info *pi) {
    bool bEOI = false;
    
    while(!pIn->isEndOfStream()) {
      /* Get next marker */
      int nMarker = pIn->readMarker();
      
      switch(nMarker) {
        case 0xDA: /* Start Of Scan */
        {
          Scan Scan;

          _ParseSOS(pIn,pRGB,pi,&Scan);
          
          _DecodeScan(pIn,pRGB,pi,&Scan);
          break;
        }
        case 0xD9: /* End Of Image */
        {
          /* Under normal circumstances we'll never get here, since the 
             EOI is processed by _DecodeScan(). Only if there's no scan (unlikely)
             we'll get here */
          bEOI = true;
          break;
        }
        default:
        {
          _InterpretMarker(pIn,pRGB,pi,nMarker);
          break;
        }
      }
      
      /* Reached end of image? */
      if(bEOI) break;
    }
  }
  
  void Decoder::_DecodeScan(InputStream *pIn,RGBImage *pRGB,Info *pi,Scan *ps) {
    /* Calculate number of restart intervals */
    int nRestartIntervalsLeft = pi->nTotalMCUs / pi->nRestartInterval;
    if(pi->nTotalMCUs % pi->nRestartInterval) nRestartIntervalsLeft++;
    
    /* Allocate memory for MCU buffer */
    if(pi->pcMCUBuffer != NULL) delete [] pi->pcMCUBuffer;
    pi->pcMCUBuffer = new uint8[g_MCULayouts[pi->nMCULayout].nMCUXSize * 
                                g_MCULayouts[pi->nMCULayout].nMCUYSize * 3]; 
    
    /* Continue until out of restart intervals */
    while(nRestartIntervalsLeft > 0) {
      /* Decode restart interval */
      int nMarker = _DecodeRestartInterval(pIn,pRGB,pi,ps);
      
      if(nMarker == 0xD9) {
        /* End Of Image marker */
        break; 
      }
      else if(nMarker >= 0xD0 && nMarker <= 0xD7) {
        /* Restart marker. If we were pedantic we'd check that they come in 
           the right order :P */
        nRestartIntervalsLeft--; 
      }
      else if(nMarker == 0) {
        /* End of stream */
        break;
      }
      else {
        /* Unexpected marker */
        throw DecodeError(ERR_UNEXPECTED_MARKER);
      }
    }    
  }
  
  int Decoder::_DecodeRestartInterval(InputStream *pIn,RGBImage *pRGB,Info *pi,Scan *ps) {
    /* Reset decoder */
    for(int i=0;i<ps->nNumComponents;i++)
      ps->Components[i].nDCPred = 0;     

    /* 8x8 block buffer */
    ALIGN16 short nBlock[64];
      
    /* Start decoding bitstream */
    pIn->flushBitStream();
    
    int nMCUsLeft = pi->nRestartInterval;
    
    while(nMCUsLeft > 0) {
      /* For each data unit in the MCU... */
      for(int i=0;i<g_MCULayouts[pi->nMCULayout].nDataUnitsPerMCU;i++) {
        /* Decode a 8x8 DCT block */
        memset(nBlock,0,sizeof(nBlock));
        
        int nComponent = g_MCULayouts[pi->nMCULayout].nDUComponents[i];
      
        /* Resolve huffman tables */
        HuffmanTable *pDC = ps->Components[nComponent].pDCTable;
        HuffmanTable *pAC = ps->Components[nComponent].pACTable;
        
        /* Decode DC value */
        nBlock[0] = ps->Components[nComponent].nDCPred + pDC->decode(pIn,NULL);
        ps->Components[nComponent].nDCPred = nBlock[0];
        
        /* Decode AC values */
        for(int j=1;j<64;j++) {
          int nLeadingZeros;
          int n = pAC->decode(pIn,&nLeadingZeros);
          
          if(n == JPG_HUFFMAN_EOB) {
            /* End-of-block symbol (remaining AC values are zero) */
            break;          
          }
          else if(n == JPG_HUFFMAN_ZRL) {
            /* 16 zeros */
            j += 15; /* last one is added implicitly by the for-loop */
          }
          else {
            /* 0-15 zeros followed by a value */
            j += nLeadingZeros;
            
            /* AC coefficients are stored in diagonal zig-zag order, starting at top left */
            nBlock[m_nZZOrder[j]] = n;            
          }
        }
        
        /* Resolve quantization matrix */
        int nQTableIndex = pi->Components[nComponent].nQuantizationTable;
        
        QuantizationTable *pQ = &pi->QuantizationTables[nQTableIndex];
       
        /* Perform dequantization -- note that we already de-zigzagged 
           the quantization matrix */        
        #if defined(EMMINTRIN)
          /* Use SIMD (only slightly faster since a lot of time is wasted loading registers) */
          const __m128i *pQRow = (const __m128i *)pQ->nTable;
          __m128i *pBRow = (__m128i *)nBlock;          
          
          for(int j=0;j<8;j++)
            pBRow[j] = _mm_mullo_epi16(pBRow[j],pQRow[j]);
        #else
          /* No SIMD */
          for(int j=0;j<64;j++)
            nBlock[j] *= pQ->nTable[j];
        #endif
        
        /* Perform IDCT on block */
        _InverseDCT(nBlock);        
        
        /* Merge 8x8 block into MCU */
        _UpdateMCU(pi,i,nBlock);
      }      
      
      /* Merge MCU into image */
      _MergeMCU(pRGB,pi);
      
      /* Next MCU... */
      nMCUsLeft--;      
      pi->nMCUIndex++;
    }    

    /* Go to next marker and return its value */
    if(pIn->getEmbeddedMarker() != 0) return pIn->getEmbeddedMarker();    
    return pIn->readMarker();    
  }
  
  void Decoder::_ParseDRI(InputStream *pIn,RGBImage *pRGB,Info *pi) {
    /* The DRI marker defines how many MCUs should be decoded at a 
       time in a single run (the restart interval). This marker is 
       optional, if it's not present we'll just do all MCUs in one run */       
    int nLen = pIn->readWord();     
    pi->nRestartInterval = pIn->readWord();
  }
  
  void Decoder::_ParseSOF(InputStream *pIn,RGBImage *pRGB,Info *pi) {
    int nLen = pIn->readWord(); 
    
    int nBitsPerComponent = pIn->readByte();    
    if(nBitsPerComponent != 8) 
      throw DecodeError(ERR_MUST_BE_8BIT);
      
    /* Get image dimensions */
    pRGB->nHeight = pIn->readWord();
    pRGB->nWidth = pIn->readWord();    
    if(pRGB->nHeight == 0 || pRGB->nWidth == 0)
      throw DecodeError(ERR_INVALID_FRAME);
    
    int nNumComponents = pIn->readByte();
    
    /* For each component... */
    for(int i=0;i<nNumComponents;i++) {
      /* Get component index */
      int nId = pIn->readByte() - 1;
      if(nId < 0 || nId >= 3) 
        throw DecodeError(ERR_INVALID_FRAME);
           
      /* Get component info */      
      int nSamplingFactors = pIn->readByte();
      pi->Components[nId].nXSampleFactor = (nSamplingFactors & 0xf0) >> 4;      
      pi->Components[nId].nYSampleFactor = nSamplingFactors & 0x0f;
      pi->Components[nId].nQuantizationTable = pIn->readByte();      
    }
    
    /* Figure out the MCU layout... */
    pi->nMCULayout = _LookupMCULayout(pi->Components,nNumComponents);

    /* Calculate total number of MCUs and default restart interval */
    pi->nHorizontalMCUs = pRGB->nWidth / g_MCULayouts[pi->nMCULayout].nMCUXSize;
    pi->nVerticalMCUs = pRGB->nHeight / g_MCULayouts[pi->nMCULayout].nMCUYSize;
    
    if(pRGB->nWidth % g_MCULayouts[pi->nMCULayout].nMCUXSize) pi->nHorizontalMCUs++;
    if(pRGB->nHeight % g_MCULayouts[pi->nMCULayout].nMCUYSize) pi->nVerticalMCUs++;
    
    pi->nTotalMCUs = pi->nHorizontalMCUs * pi->nVerticalMCUs;
    pi->nRestartInterval = pi->nTotalMCUs;  
  }
  
  void Decoder::_ParseDHT(InputStream *pIn,RGBImage *pRGB,Info *pi) {
    int nLen = pIn->readWord() - 2; 

    /* Read huffman tables */
    while(nLen > 0) {
      int n = pIn->readByte();
      nLen--;

      int nId = n & 0x0f;
      if(nId < 0 || nId >= 3)
        throw DecodeError(ERR_INVALID_HUFFMAN_TABLE_ID);
            
      HuffmanTable *pTable;
            
      if((n & 0xf0) >> 4)
        pTable = &pi->ACHuffmanTables[nId];
      else
        pTable = &pi->DCHuffmanTables[nId];
      
      for(int i=0;i<16;i++)
        pTable->nL[i] = pIn->readByte();
      nLen -= 16;
      
      for(int i=0;i<16;i++) {
        if(pTable->nL[i] > 0) pTable->nMaxBits = i + 1;
      
        for(int j=0;j<pTable->nL[i];j++) {
          pTable->nTable[i*256+j] = pIn->readByte();
        }
        nLen -= pTable->nL[i];
      }

      /* Build tree */
      pTable->buildTree();
        
      pTable->bDefined = true;
    }
  }
  
  void Decoder::_ParseDQT(InputStream *pIn,RGBImage *pRGB,Info *pi) {
    int nLen = pIn->readWord() - 2; 
    
    /* Read quantization table definitions */
    while(nLen > 0) {
      int n = pIn->readByte();
      nLen--;
      
      int nId = n & 0x0f;
      if(nId < 0 || nId >= 3)
        throw DecodeError(ERR_INVALID_QUANTIZATION_TABLE);
      
      int nPrecission = (n & 0xf0) >> 4;      
      if(nPrecission != 0)
        throw DecodeError(ERR_INVALID_QUANTIZATION_TABLE);
              
      pi->QuantizationTables[nId].bDefined = true;
      
      for(int i=0;i<64;i++) {
        /* Un-zigzag the order */
        pi->QuantizationTables[nId].nTable[m_nZZOrder[i]] = pIn->readByte();
      }
      
      nLen -= 64;
    }
  }  
  
  void Decoder::_ParseSOS(InputStream *pIn,RGBImage *pRGB,Info *pi,Scan *ps) {
    int nLen = pIn->readWord();
    
    ps->nNumComponents = pIn->readByte();
    
    /* For each component... */
    for(int i=0;i<ps->nNumComponents;i++) {
      int nId = pIn->readByte() - 1;
      if(nId < 0 || nId >= 3)
        throw DecodeError(ERR_INVALID_SCAN);
      
      ps->Components[i].nIndex = nId;
      
      int nEncodingTables = pIn->readByte();
      
      ps->Components[i].nACTable = nEncodingTables & 0x0f;
      ps->Components[i].nDCTable = (nEncodingTables & 0xf0) >> 4;
      
      /* Resolve huffman tables */
      ps->Components[i].pACTable = &pi->ACHuffmanTables[ps->Components[i].nACTable];
      ps->Components[i].pDCTable = &pi->DCHuffmanTables[ps->Components[i].nDCTable];
      
      if(!ps->Components[i].pACTable->bDefined ||
         !ps->Components[i].pDCTable->bDefined)
        throw DecodeError(ERR_INVALID_SCAN);         
    }        
    
    ps->nDCTStart = pIn->readByte();
    ps->nDCTEnd = pIn->readByte();
    
    int n = pIn->readByte();
    
    ps->nAH = n & 0x0f;
    ps->nAL = (n & 0xf0) >> 4;
  }
  
  void Decoder::_MergeMCU(RGBImage *pRGB,Info *pi) {
    /* Merge MCU into the image while making sure it's clipped against the right/bottom edges */
    int nMCUClippedXSize = g_MCULayouts[pi->nMCULayout].nMCUXSize;
    int nMCUClippedYSize = g_MCULayouts[pi->nMCULayout].nMCUYSize;;
    
    int nX = (pi->nMCUIndex % pi->nHorizontalMCUs) * nMCUClippedXSize;
    int nY = (pi->nMCUIndex / pi->nHorizontalMCUs) * nMCUClippedYSize;
    
    if(nX + g_MCULayouts[pi->nMCULayout].nMCUXSize > pRGB->nWidth) 
      nMCUClippedXSize -= (nX + g_MCULayouts[pi->nMCULayout].nMCUXSize) - pRGB->nWidth;
      
    if(nY + g_MCULayouts[pi->nMCULayout].nMCUYSize > pRGB->nHeight)
      nMCUClippedYSize -= (nY + g_MCULayouts[pi->nMCULayout].nMCUYSize) - pRGB->nHeight;

    int nPitchDest = pRGB->nWidth * 3;    
    int nPitchSrc = g_MCULayouts[pi->nMCULayout].nMCUXSize * 3;
    uint8 *pcDest = &pRGB->pc[nPitchDest * nY + nX * 3];
    uint8 *pcSrc = pi->pcMCUBuffer;
    
    for(int i=0;i<nMCUClippedYSize;i++) {
      uint8 *pcD = pcDest;
      uint8 *pcS = pcSrc;
      
      for(int j=0;j<nMCUClippedXSize;j++) {
        /* Convert from YCBCR to RGB */
        int r = pcS[0] + m_nRTable[pcS[2]];
        int g = pcS[0] + m_nGTable[((int)pcS[1]<<8)+pcS[2]];
        int b = pcS[0] + m_nBTable[pcS[1]];    
          
        pcD[0] = BYTECLAMP(r);
        pcD[1] = BYTECLAMP(g);
        pcD[2] = BYTECLAMP(b);          
      
        pcD+=3;
        pcS+=3;
      }
      pcDest += nPitchDest;
      pcSrc += nPitchSrc;
    }        
  }
  
  void Decoder::_UpdateMCU8x8(uint8 *pcMCUBuffer,int nMCUXSize,int nMCUYSize,short *pnBlock,int nComponent,int nXOffset,int nYOffset) {
    /* This function merges an unscaled 8x8 block into an arbitrary sized MCU/channel */
    uint8 *pc = &pcMCUBuffer[(nYOffset * nMCUXSize + nXOffset) * 3 + nComponent];
    short *pnD = pnBlock;
    
    for(int i=0;i<8;i++) {
      uint8 *pcC = pc;
      for(int j=0;j<8;j++) {
        int n = BYTECLAMP(*pnD);
        pnD++;
        
        pcC[0] = n;
        
        pcC += 3;
      }
      pc += nMCUXSize * 3;
    }
  }

  void Decoder::_UpdateMCU16x8(uint8 *pcMCUBuffer,short *pnBlock,int nComponent) {
    /* This function stretches a 8x8 block to 16x8 */
    uint8 *pc = &pcMCUBuffer[nComponent];
    short *pnD = pnBlock;
    
    for(int i=0;i<8;i++) {
      uint8 *pcC = pc;
      for(int j=0;j<8;j++) {
        int n = BYTECLAMP(*pnD);
        pnD++;
        
        pcC[0] = n; pcC[3] = n;
        
        pcC += 6;
      }
      pc += 48;
    }
  }

  void Decoder::_UpdateMCU16x16(uint8 *pcMCUBuffer,short *pnBlock,int nComponent) {
    /* This function stretches a 8x8 block to 16x16 */
    uint8 *pc = &pcMCUBuffer[nComponent];
    short *pnD = pnBlock;
    
    for(int i=0;i<8;i++) {
      uint8 *pcC = pc;
      for(int j=0;j<8;j++) {
        int n = BYTECLAMP(*pnD);
        pnD++;
                
        pcC[0] = n; pcC[3] = n; pcC[48] = n; pcC[51] = n;
        
        pcC += 6;
      }
      pc += 96;
    }
  }

  void Decoder::_UpdateMCU(Info *pi,int nDataUnit,short *pnBlock) {
    /* Updates the current MCU with the given data unit */
    switch(pi->nMCULayout) {
      case JPG_MCU_LAYOUT_YCBCR_420:
        switch(nDataUnit) {
          case 0: _UpdateMCU8x8(pi->pcMCUBuffer,16,16,pnBlock,0,0,0); break;
          case 1: _UpdateMCU8x8(pi->pcMCUBuffer,16,16,pnBlock,0,8,0); break;
          case 2: _UpdateMCU8x8(pi->pcMCUBuffer,16,16,pnBlock,0,0,8); break;
          case 3: _UpdateMCU8x8(pi->pcMCUBuffer,16,16,pnBlock,0,8,8); break;
          case 4: _UpdateMCU16x16(pi->pcMCUBuffer,pnBlock,1); break;
          case 5: _UpdateMCU16x16(pi->pcMCUBuffer,pnBlock,2); break;
        }
        break;
      case JPG_MCU_LAYOUT_YCBCR_422:
        switch(nDataUnit) {
          case 0: _UpdateMCU8x8(pi->pcMCUBuffer,16,8,pnBlock,0,0,0); break;
          case 1: _UpdateMCU8x8(pi->pcMCUBuffer,16,8,pnBlock,0,8,0); break;
          case 2: _UpdateMCU16x8(pi->pcMCUBuffer,pnBlock,1); break;
          case 3: _UpdateMCU16x8(pi->pcMCUBuffer,pnBlock,2); break;
        }
        break;
      case JPG_MCU_LAYOUT_YCBCR_444:
        switch(nDataUnit) {
          case 0: _UpdateMCU8x8(pi->pcMCUBuffer,8,8,pnBlock,0,0,0); break;
          case 1: _UpdateMCU8x8(pi->pcMCUBuffer,8,8,pnBlock,1,0,0); break;
          case 2: _UpdateMCU8x8(pi->pcMCUBuffer,8,8,pnBlock,2,0,0); break;
        }
        break;
    }
  }
  
  void Decoder::_InitRGBTables(void) {            
    /* Init tables for fast YCBCR-to-RGB conversion */
    for(int i=0;i<256;i++) {
      m_nRTable[i] = ((91881 * (i - 128)) >> 16);  /* Min=-180, Max=180 */
      m_nBTable[i] = ((116129 * (i - 128)) >> 16); /* Min=-227, Max=227 */
      
      for(int j=0;j<256;j++)
        m_nGTable[(i<<8) + j] = ((-22553 * (i - 128) - 46801 * (j - 128)) >> 16); /* Min=-135, Max=135 */
    }
  }
  
  void Decoder::_InitZZOrdering(void) {
    /* DCT matrices are stored in a diagonal zig-zag order, starting
       at the upper left corner (the DC value). This function initializes
       a table which convert a sequential index to a zig-zag'd index */
    int q = 0;
    
    for(int x=0;x<8;x++) {
      for(int i=0;i<8;i++) {
        if(x % 2 == 0) {
          int xx = i;
          int yy = x - i;
          if(yy < 0) break;          
          
          m_nZZOrder[q++] = xx + yy * 8;
        }
        else {
          int xx = x - i;
          int yy = i;
          if(xx < 0) break;          
          
          m_nZZOrder[q++] = xx + yy * 8;
        }
      }
    }
    for(int y=1;y<8;y++) {
      for(int i=0;i<8;i++) {
        if(y % 2 == 0) {
          int xx = 7 - i;
          int yy = y + i;
          if(yy >= 8) break;

          m_nZZOrder[q++] = xx + yy * 8;
        }
        else {
          int xx = y + i;
          int yy = 7 - i;
          if(xx >= 8) break;

          m_nZZOrder[q++] = xx + yy * 8;
        }
      }
    }
  }
  
  /*===========================================================================
  2D IDCTs are calculated in two passes of 1D IDCTs, first on the rows and 
  then the columns.
  ===========================================================================*/
  template <int _Half,int _FinalRShift,int _FinalOffset>
  void _IDCT_1D(const short *pnIn,short *pnOut) {
    /* 1D IDCT butterfly solution (LLM) */
    int p, n;
    
    int x0 = pnIn[0] << 9;
    int x1 = pnIn[1] << 7;
    int x2 = pnIn[2];
    int x3 = pnIn[3] * 181;
    int x4 = pnIn[4] << 9;
    int x5 = pnIn[5] * 181;
    int x6 = pnIn[6];
    int x7 = pnIn[7] << 7;
    
    n = 277*(x6+x2); p = x6;
    x6 = n + 392*x2;
    x2 = n - 946*p;
    
    p = x0 + x4; n = x0 - x4;
    x0 = p + x6 + _Half;
    x4 = n + x2 + _Half;
    x6 = p - x6 + _Half;
    x2 = n - x2 + _Half;
    
    p = x1 + x7; n = x1 - x7;
    x1 = p + x3;
    x7 = n + x5;
    x3 = p - x3;
    x5 = n - x5;
    
    n = 251*(x5+x3); p = x5;
    x5 = (n - 201*x3) >> 6;
    x3 = (n - 301*p) >> 6;
    
    n = 213*(x1+x7); p = x1;
    x1 = (n - 71*x7) >> 6;
    x7 = (n - 355*p) >> 6;
    
    pnOut[0] = (short)((x0 + x1) >> _FinalRShift) + _FinalOffset;
    pnOut[8] = (short)((x4 + x5) >> _FinalRShift) + _FinalOffset;
    pnOut[16] = (short)((x2 + x3) >> _FinalRShift) + _FinalOffset;
    pnOut[24] = (short)((x6 + x7) >> _FinalRShift) + _FinalOffset;
    pnOut[32] = (short)((x6 - x7) >> _FinalRShift) + _FinalOffset;
    pnOut[40] = (short)((x2 - x3) >> _FinalRShift) + _FinalOffset;
    pnOut[48] = (short)((x4 - x5) >> _FinalRShift) + _FinalOffset;
    pnOut[56] = (short)((x0 - x1) >> _FinalRShift) + _FinalOffset;
  }

  void Decoder::_InverseDCT(short *pnBlock) {  
    /* Calculate inverse DCT */
    short nTemp[64];

    /* First perform 1D IDCT on each row */
    for(int i=0;i<8;i++)
      _IDCT_1D<256,9,0>(&pnBlock[i*8],&nTemp[i]);
    
    /* Now do the same 1D IDCT on each column -- first pass generated a transposed
       matrix so we can use the same algorithm again */  
    for(int i=0;i<8;i++)
      _IDCT_1D<2048,12,128>(&nTemp[i*8],&pnBlock[i]);
  }  
    
}
