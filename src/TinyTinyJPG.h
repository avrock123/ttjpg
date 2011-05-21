/*
TinyTinyJPG library
Copyright (c) 2011, Rasmus Neckelmann (neckelmann@gmail.com)
All rights reserved.

This source code is released under the terms of the BSD license; see the 
file LICENSE.TXT for details.
*/
#ifndef __TINYTINYJPG_H__
#define __TINYTINYJPG_H__

#include <stdio.h>
#include <string.h>

/*=============================================================================
Definitions
=============================================================================*/
#if defined(VLD) && defined(_DEBUG)
  #include <vld.h>
#endif

#if defined(_MSC_VER) 
  #define ALIGN16 __declspec(align(16))
#else
  #define ALIGN16
#endif

/* Special huffman symbol codes */
#define JPG_HUFFMAN_EOB             0
#define JPG_HUFFMAN_ZRL             -1000000

/* Known and supported MCU layouts */
#define JPG_MCU_LAYOUT_YCBCR_444    0
#define JPG_MCU_LAYOUT_YCBCR_422    1
#define JPG_MCU_LAYOUT_YCBCR_420    2
#define JPG_NUM_MCU_LAYOUTS         3

namespace ttjpg {

  /* Forward declarations */
  class InputStream;

  /*===========================================================================
  Types
  ===========================================================================*/
  /* Unsigned types */
  typedef unsigned char uint8;
  typedef unsigned short uint16;
  typedef unsigned long uint32;
  
  /* Error codes */
  enum ErrorCode {
    ERR_NONE,ERR_UNSUPPORTED_MCU_LAYOUT,ERR_INVALID_HUFFMAN_TREE,
    ERR_INVALID_HUFFMAN_CODE,ERR_NOT_JPEG,ERR_NO_FRAME,ERR_UNEXPECTED_MARKER,
    ERR_MUST_BE_8BIT,ERR_INVALID_FRAME,ERR_INVALID_HUFFMAN_TABLE_ID,
    ERR_INVALID_QUANTIZATION_TABLE,ERR_INVALID_SCAN,
    ERR_DYNAMIC_ARRAY_TOO_LARGE,ERR_STREAM_FAILURE
  };

  /* Simple 24-bit RGB (8 bits per component) images. Note that the pixel data
     memory is released when the object runs out of scope */
  struct RGBImage {
    RGBImage() {
      /* Tidy */
      pc = NULL;
      nWidth = nHeight = 0;
    }
    ~RGBImage() {
      /* Automatic release memory */
      release();
    }
    
    /* Data */
    uint8 *pc;
    int nWidth,nHeight;
    
    /* Methods */
    void release(void) {
      if(pc != NULL) delete [] pc;
      pc = NULL;
    }
  };
  
  /* Component information as described in the JPEG file */
  struct ComponentInfo {
    ComponentInfo() {
      /* Defaults */
      nXSampleFactor = 1;
      nYSampleFactor = 1;
      nQuantizationTable = 0;
    }
    
    int nXSampleFactor;
    int nYSampleFactor;
    int nQuantizationTable;
  };
  
  /* 8x8 quantization table (table is in zig-zag order) */
  struct QuantizationTable {
    QuantizationTable() {
      /* Defaults */
      bDefined = false;
      memset(nTable,0,sizeof(nTable));
    }
    
    bool bDefined;
    ALIGN16 short nTable[64];
  };
  
  /* Node of the huffman decoding tree */  
  struct HuffmanNode {
    HuffmanNode() {
      /* Defaults */
      pChildren[0] = pChildren[1] = NULL;
      nSymbol = -1;
    }
    
    int nSymbol;
    HuffmanNode *pChildren[2];
  };
  
  /* Dynamic array. Allocations are done in buckets. */
  template<typename _T,int _AllocSize,int _MaxAllocs>
  struct DynamicArray {
    DynamicArray() {
      _Tidy();
    }
    ~DynamicArray() {
      clear();
    }
  
    _T Static[_AllocSize];
    _T *p[_MaxAllocs];
    int nAllocs,nSize;
    
    /* Methods */
    void clear(void) {
      for(int i=0;i<nAllocs;i++)
        delete [] p[i];
      _Tidy();
    }
    
    void grow(int n) {
      int nNewSize = nSize + n;
      
      if(nNewSize > (nAllocs + 1) * _AllocSize) {
        int nNeeded = (nNewSize / _AllocSize) - 1;
        if(nNewSize % _AllocSize) nNeeded++;
      
        if(nAllocs >= _MaxAllocs)
          throw DecodeError(ERR_DYNAMIC_ARRAY_TOO_LARGE);
        
        for(int i=nAllocs;i<nNeeded;i++)
          p[i] = new _T[_AllocSize];
        
        nAllocs = nNeeded;
      }
            
      nSize = nNewSize;
    }
    
    int size(void) const {
      return nSize;
    }
    
    _T & operator [] (int i) {
      int a = i / _AllocSize;
      if(a == 0) return Static[i % _AllocSize];
      return p[a][(i % _AllocSize) - 1];
    }
    
    /* Local functions */
    void _Tidy(void) {
      memset(p,0,sizeof(_T *) * _MaxAllocs);
      nAllocs = nSize = 0;
    }
  };
  
  /* Huffman table and associated decoding tree built from it */
  struct HuffmanTable {
    HuffmanTable() {
      /* Defaults */
      bDefined = false;
      memset(nTable,0,sizeof(nTable));
      memset(nL,0,sizeof(nL));
      
      nNextFreeNode = 0;
      nNumFreeNodes = 0;
      nMaxBits = 0;
    }
    ~HuffmanTable() {
    }
    
    bool bDefined;
    
    int nMaxBits;     /* Maximum code length (i.e. height of tree) */    
    int nL[16];       
    int nTable[4096];
    
    int nNextFreeNode;
    int nNumFreeNodes;
    DynamicArray<HuffmanNode,1024,1024> Tree;
    
    /* Methods */
    void buildTree(void);
    void buildLookupTables(HuffmanNode *p,int nBits,uint32 nCode);
    int decode(InputStream *pIn,int *pnLeadingZeros);    
  };
  
  /* JPEG information (internal per-decode context) */
  struct Info {
    Info() {
      pcMCUBuffer = NULL;
    }
    ~Info() {
      if(pcMCUBuffer != NULL) delete [] pcMCUBuffer;
    }
  
    ComponentInfo Components[3];
    QuantizationTable QuantizationTables[4];
    HuffmanTable ACHuffmanTables[4];
    HuffmanTable DCHuffmanTables[4];
    
    int nRestartInterval;
    
    int nMCULayout;
    int nMCUIndex;
    int nTotalMCUs;
    
    int nHorizontalMCUs,nVerticalMCUs;
    
    uint8 *pcMCUBuffer;
  };
  
  /* Component of a scan */
  struct ScanComponent {
    HuffmanTable *pACTable; /* Huffman decoding of AC values */
    HuffmanTable *pDCTable; /* Huffman decoding of DC values */   
    int nDCPred;            /* DC predictor */
  };  
  
  /* Internal per-decode scan context */
  struct Scan {        
    ScanComponent Components[3];     
  };    

  /*===========================================================================
  Exception class
  ===========================================================================*/
  class DecodeError {
    public:
      DecodeError(ErrorCode _ErrorCode) : m_ErrorCode(_ErrorCode) {}
      ~DecodeError() {}
      
      /* Methods */
      ErrorCode getError(void) {return m_ErrorCode;}
    
    private:
      /* Data */
      ErrorCode m_ErrorCode;
  };

  /*===========================================================================
  InputStream abstract class interface
  
    Applications should implement this class in order to stream data to the
    decoder.
  ===========================================================================*/
  class InputStream {
    public:
      InputStream() {
        flushBitStream();
      }
      virtual ~InputStream() {
      }
      
      /*-----------------------------------------------------------------------
      Virtual interface 
      -----------------------------------------------------------------------*/     
      /* isEndOfStream() should return true if end-of-stream has been hit 
         crossed */
      virtual bool isEndOfStream(void) = 0;
      
      /* readByte() should read a single byte from the stream and return it */
      virtual int readByte(void) = 0;
     
      /* skip() should skip the next nBytes bytes in stream */
      virtual void skip(int nBytes) = 0;
      
      /*-----------------------------------------------------------------------
      Methods
      -----------------------------------------------------------------------*/     
      /* Common stream helper functions */
      int readWord(void);
      int readMarker(void);     
      uint32 readBits(int nNum);
      void flushBitStream(void);
      
      int getEmbeddedMarker(void) {return m_nEmbeddedMarker;}
      
    private:
      /*-----------------------------------------------------------------------
      Data
      -----------------------------------------------------------------------*/          
      uint32 m_nBitBuf;
      int m_nBits;      
      int m_nEmbeddedMarker;
  };  
  
  /*===========================================================================
  StdioInputStream class
  
    stdio (fopen(), etc) implementation of the InputStream interface
  ===========================================================================*/
  class StdioInputStream : public InputStream {
    public:
      StdioInputStream(FILE *fp,bool bCloseFile = false) 
        : m_fp(fp),
          m_bCloseFile(bCloseFile),
          m_nBufferedBytes(0),
          m_nNextByte(0),
          m_bEOS(false) {
        /* Use supplied stdio stream */
      }
      StdioInputStream(const char *pcFileName)
        : m_bCloseFile(true),
          m_nBufferedBytes(0),
          m_nNextByte(0),
          m_bEOS(false) {
        /* Open file */
        m_fp = fopen(pcFileName,"rb");
      }      
      virtual ~StdioInputStream() {
        /* Close handle? */
        if(m_bCloseFile && m_fp != NULL)
          fclose(m_fp);
      }
      
      /* InputStream implementation */
      virtual bool isEndOfStream(void) {
        return m_bEOS;
      }
      
      virtual int readByte(void) {
        if(m_nBufferedBytes == 0) {        
          /* Need to fill the buffer */
          if(m_fp == NULL) throw DecodeError(ERR_STREAM_FAILURE);
          
          m_nBufferedBytes = (int)fread(m_cBuffer,1,sizeof(m_cBuffer),m_fp);
          
          if(ferror(m_fp) != 0) throw DecodeError(ERR_STREAM_FAILURE);
          
          if(m_nBufferedBytes == 0) {
            m_bEOS = true;
            return 0;
          }
          
          m_nNextByte = 0;
        }
        
        /* Grab next byte and return it */
        int n = m_cBuffer[m_nNextByte];
        m_nNextByte++;
        m_nBufferedBytes--;
        return n;
      }
      
      virtual void skip(int nBytes) {
        /* Skip next bytes */
        if(m_nBufferedBytes >= nBytes) {
          m_nBufferedBytes -= nBytes;
          m_nNextByte += nBytes;
        }
        else {
          int nRemaining = nBytes - m_nBufferedBytes;
        
          m_nBufferedBytes = 0;
          m_nNextByte = 0;
          
          fseek(m_fp,nRemaining,SEEK_CUR);

          if(ferror(m_fp) != 0) throw DecodeError(ERR_STREAM_FAILURE);
        }
      }      
      
    private:
      /* Data */
      FILE *m_fp;
      bool m_bCloseFile;
      bool m_bEOS;
      
      /* Buffer */
      int m_nBufferedBytes;
      int m_nNextByte;
      uint8 m_cBuffer[16384];
  };
  
  /*===========================================================================
  MemoryInputStream class
  
    Implementation of the InputStream interface for reading JPEG files directly
    from a memory buffer.
  ===========================================================================*/
  class MemoryInputStream : public InputStream {
    public:
      MemoryInputStream(const uint8 *pc,int nTotalBytes) 
        : m_pc(pc),
          m_nTotalBytes(nTotalBytes),
          m_nNextByte(0),
          m_bEOS(false) {
      }
      virtual ~MemoryInputStream() {
      }
      
      /* InputStream implementation */
      virtual bool isEndOfStream(void) {
        return m_bEOS;
      }
      
      virtual int readByte(void) {
        if(m_nNextByte == m_nTotalBytes) {        
          m_bEOS = true;
          return 0;
        }
        
        /* Grab next byte and return it */
        return m_pc[m_nNextByte++];
      }
      
      virtual void skip(int nBytes) {
        /* Skip next bytes */
        m_nNextByte += nBytes;
        
        if(m_nNextByte > m_nTotalBytes) {
          m_bEOS = true;
          m_nNextByte = m_nTotalBytes;
        }
      }      
      
    private:
      /* Data */
      bool m_bEOS;
      int m_nTotalBytes;
      int m_nNextByte;
      const uint8 *m_pc;
  };
  
  /*===========================================================================
  Decoder class
  
    This class defines the decoding functionality of the library.
  ===========================================================================*/
  class Decoder {
    public:
      Decoder();
      ~Decoder();
      
      /*-----------------------------------------------------------------------
      Methods
      -----------------------------------------------------------------------*/     
      void readImageInfo(InputStream *pIn,RGBImage *pRGB);
      void readImage(InputStream *pIn,RGBImage *pRGB);    
      void decodeScan(InputStream *pIn,RGBImage *pRGB,Info *pi,Scan *ps);
      
    private:
      /*-----------------------------------------------------------------------
      Data
      -----------------------------------------------------------------------*/     
      int m_nZZOrder[64];
      
      short m_nRTable[256];
      short m_nGTable[65536];
      short m_nBTable[256];      
    
      /*-----------------------------------------------------------------------
      Local functions
      -----------------------------------------------------------------------*/     
      void _InterpretMarker(InputStream *pIn,RGBImage *pRGB,Info *pi,int nMarker);
      void _DecodeFrame(InputStream *pIn,RGBImage *pRGB,Info *pi);
      int _DecodeRestartInterval(InputStream *pIn,RGBImage *pRGB,Info *pi,Scan *ps);
      
      void _Parse(InputStream *pIn,RGBImage *pRGB,Info *pi,bool bGetInfoOnly);
      void _ParseSOF(InputStream *pIn,RGBImage *pRGB,Info *pi);
      void _ParseDHT(InputStream *pIn,RGBImage *pRGB,Info *pi);
      void _ParseDQT(InputStream *pIn,RGBImage *pRGB,Info *pi);
      void _ParseSOS(InputStream *pIn,RGBImage *pRGB,Info *pi,Scan *ps);
      void _ParseDRI(InputStream *pIn,RGBImage *pRGB,Info *pi);
      
      void _InverseDCT(short *pnBlock);
      
      void _MergeMCU(RGBImage *pRGB,Info *pi);      
      void _UpdateMCU(Info *pi,int nDataUnit,short *pnBlock);
      
      void _UpdateMCU8x8(uint8 *pcMCUBuffer,int nMCUXSize,int nMCUYSize,short *pnBlock,int nComponent,int nXOffset,int nYOffset);
      void _UpdateMCU16x8(uint8 *pcMCUBuffer,short *pnBlock,int nComponent);
      void _UpdateMCU16x16(uint8 *pcMCUBuffer,short *pnBlock,int nComponent);
      
      void _InitZZOrdering(void);
      void _InitRGBTables(void);            
  };

}

#endif
