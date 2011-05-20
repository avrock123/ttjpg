/*
TinyTinyJPG library
Copyright (c) 2011, Rasmus Neckelmann (neckelmann@gmail.com)
All rights reserved.

This source code is released under the terms of the BSD license; see the 
file LICENSE.TXT for details.
*/

/*

  IMPORTANT: 
  
        (1)  In order for this to run succesfully you'll need the proper
             files in ./testdata/
             
             These files are provided as a seperate download.       
             
        (2)  unit++ required to compile.

*/

#include "unitpp/unit++.h"

#include "../src/TinyTinyJPG.h"

#define TTJPG_BENCHMARK
#define TTJPG_CHECK_PIXELS
//#define JPGD

#if defined(JPGD)
  /* Just drag it in like a whale */
  #include "jpgd.cpp"
#endif

#if defined(TTJPG_BENCHMARK)
  #if defined(WIN32)
    #include <windows.h>
  #else
    #error "benchmark only on win32"
  #endif
#endif

namespace tests {

  using namespace unitpp;
  
  #if defined(TTJPG_BENCHMARK)
  
  #if defined(WIN32)
    __int64 g_nPerfFreq = 0;
  #endif
  
  class Timer {
    public:
      Timer() {
        #if defined(WIN32)
          if(g_nPerfFreq == 0) {
            QueryPerformanceFrequency((LARGE_INTEGER *)&g_nPerfFreq);
            
            /* Make sure we don't switch CPUs (can mess up timing) */
            SetThreadAffinityMask(GetCurrentThread(),1);
          }
        #endif
      }
      ~Timer() {
      } 
      
      void start(void) {
        #if defined(WIN32)
          QueryPerformanceCounter((LARGE_INTEGER *)&m_nStart);
        #endif
      }
      
      int stop(void) {
        /* Return elapsed time in milliseconds */
        #if defined(WIN32)
          __int64 n;
          QueryPerformanceCounter((LARGE_INTEGER *)&n);
          return (int)(((double)(n - m_nStart) / (double)g_nPerfFreq) * 1000.0);
        #endif
      }
      
      #if defined(WIN32)
        __int64 m_nStart;
      #endif
  };
  #endif

  class StdioInput : public ttjpg::InputStream {
    public:
      StdioInput(FILE *fp) : m_fp(fp) {}
      virtual ~StdioInput() {}
      
      virtual bool isEndOfStream(void) {
        return feof(m_fp) != 0;
      }
      
      virtual int readByte(void) {
        return fgetc(m_fp);
      }
      
      virtual void skip(int nBytes) {
        fseek(m_fp,nBytes,SEEK_CUR);
      }
        
    private:
      FILE *m_fp;
  };
  
  class BufferedStdioInput : public ttjpg::InputStream {
    public:
      BufferedStdioInput(FILE *fp) {
        fseek(fp,0,SEEK_END);
        m_nSize = ftell(fp);
        fseek(fp,0,SEEK_SET);
        m_nPos = 0;
        m_bEOS = false;
        
        m_pc = new unsigned char[m_nSize];
        
        fread(m_pc,m_nSize,1,fp);
      }
      virtual ~BufferedStdioInput() {
        delete [] m_pc;
      }
      
      virtual bool isEndOfStream(void) {
        return m_bEOS;
      }
      
      virtual int readByte(void) {
        if(m_nPos >= m_nSize) {
          m_bEOS = true;
          return 0;
        }
        
        return m_pc[m_nPos++];
      }
      
      virtual void skip(int nBytes) {
        m_nPos += nBytes;
        
        if(m_nPos > m_nSize) m_bEOS = true;
      }    
      
    private:
      unsigned char *m_pc;
      int m_nPos;
      int m_nSize;
      bool m_bEOS;
  };

  class ttjpgTest : public suite {
    public:
      ttjpgTest() : suite("ttjpgTest") {        
        main().add("test_decoder",testcase(this,"decoder",&ttjpgTest::test_decoder));
      }
      
      /* Tests */
      void test_decoder(void) {
        /* Create decoder object */
        ttjpg::Decoder JPG;               
        
        /* Test files */
        //_TestJPG(&JPG,"testdata/img_3156x3156_444.jpg",3156,3156);
        _TestJPG(&JPG,"testdata/lol.jpg",3072,2304);
        
        _TestJPG(&JPG,"testdata/img_320x240_420_restart5.jpg",320,240);
        _TestJPG(&JPG,"testdata/img_320x240_420_restart5B.jpg",320,240);
        
        _TestJPG(&JPG,"testdata/img_320x240_420.jpg",320,240);
        _TestJPG(&JPG,"testdata/img_320x240_422.jpg",320,240);
        _TestJPG(&JPG,"testdata/img_320x240_444.jpg",320,240);
        
        _TestJPG(&JPG,"testdata/img_321x241_420.jpg",321,241);
        _TestJPG(&JPG,"testdata/img_321x241_422.jpg",321,241);
        _TestJPG(&JPG,"testdata/img_321x241_444.jpg",321,241);

        #if defined(JPGD)
          //_TestJPG_jpgd("testdata/lol.jpg");
          //_TestJPG_jpgd("testdata/img_3156x3156_444.jpg");
          _TestJPG_jpgd("testdata/img_320x240_420.jpg");
        #endif
      }
      
    private:
      /* Helper functions */
      #if defined(JPGD)
        void _TestJPG_jpgd(const std::string &Name) {
          #if defined(TTJPG_BENCHMARK)
            printf("[%s] (jpgd)\n",Name.c_str());
          
            Timer Stopwatch;          
            Stopwatch.start();
                    
            int nWidth,nHeight,nActualComps;
            unsigned char *pc = jpgd::decompress_jpeg_image_from_file(Name.c_str(),&nWidth,&nHeight,&nActualComps,3);

            int nTime = Stopwatch.stop();
            printf("  %d ms\n",nTime);
          #endif
        }
      #endif
      
      void _TestJPG(ttjpg::Decoder *pDecoder,const std::string &Name,int nExpectWidth,int nExpectHeight) {
        #if defined(TTJPG_BENCHMARK)
          printf("[%s]\n",Name.c_str());
        #endif
      
        /* Just get image info at first */
        ttjpg::RGBImage Info;
        _GetJPGInfo(pDecoder,Name,&Info);
        
        assert_eq(Name + " info : width",nExpectWidth,Info.nWidth);
        assert_eq(Name + " info : height",nExpectHeight,Info.nHeight);
        
        /* Decode full image */
        ttjpg::RGBImage Img;
        _LoadJPG(pDecoder,Name,&Img);

        assert_eq(Name + " load : width",nExpectWidth,Img.nWidth);
        assert_eq(Name + " load : height",nExpectHeight,Img.nHeight);
        assert_true(Name + " load : data",Img.pc != NULL);
        
        /* Output .raw and compare it with the corrent result */
        FILE *fp = fopen((Name + ".raw").c_str(),"wb");        
        assert_true(Name + ".raw: fopen",fp != NULL);
        
        #if defined(TTJPG_CHECK_PIXELS)
          FILE *fp2 = fopen((Name + "_check.raw").c_str(),"rb");
          assert_true(Name + "_check.raw: fopen",fp2 != NULL);
        #endif        
        
        for(int i=0;i<Img.nWidth*Img.nHeight*3;i++) {
          fputc(Img.pc[i],fp);          
          
          #if defined(TTJPG_CHECK_PIXELS)
            assert_true(Name + " load: check pixels",abs(fgetc(fp2)-Img.pc[i]) < 10);
          #endif
        }
        
        fclose(fp);        
        
        #if defined(TTJPG_CHECK_PIXELS)
          fclose(fp2);        
        #endif
      }
      
      void _LoadJPG(ttjpg::Decoder *pDecoder,const std::string &Name,ttjpg::RGBImage *pImage) {
        /* Open file for input */
        FILE *fp = fopen(Name.c_str(),"rb");
        assert_true("fopen(): " + Name,fp != NULL);
        
        /* Setup stream */
        BufferedStdioInput In(fp);
        
        /* Decode JPG */        
        #if defined(TTJPG_BENCHMARK)        
          Timer Stopwatch;          
          Stopwatch.start();
        #endif
        
        pDecoder->readImage(&In,pImage);

        #if defined(TTJPG_BENCHMARK)        
          int nTime = Stopwatch.stop();
          printf("  %d ms\n",nTime);
        #endif
        
        /* Clean up */
        fclose(fp);
      }

      void _GetJPGInfo(ttjpg::Decoder *pDecoder,const std::string &Name,ttjpg::RGBImage *pImage) {
        /* Open file for input */
        FILE *fp = fopen(Name.c_str(),"rb");
        assert_true("fopen(): " + Name,fp != NULL);
        
        /* Setup stream */
        StdioInput In(fp);
        
        /* Decode JPG */
        pDecoder->readImageInfo(&In,pImage);
        
        /* Clean up */
        fclose(fp);
      }
  };
  
  /* Create test instance */
  ttjpgTest *g_pttjpgTest = new ttjpgTest;

}
