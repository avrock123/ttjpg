TinyTinyJPG v1.0
(C) 2011 Rasmus Neckelmann (neckelmann@gmail.com)
-------------------------------------------------

What is this?
-------------
  TinyTinyJPG (or ttjpg for short) is a simple C++ library for loading images 
  stored in the JPEG format. The general intention is to make it as small as
  possible by only supporting the most commonly used features of JPEG.
  
  The implementation is fairly straightforward, so this library should also
  be pretty useful if you just want to learn how JPEG works.
  
  A few optional SSE2 optimizations included (Visual C++ only), although the
  code could easily be made much faster at the cost of getting more convoluted.

Limitations
-----------
  * YCbCr-color, interleaved (single scan), baseline JPEGs only. 
  
  * Only the 3 most common chroma-subsampling modes (MCU layouts) are 
    supported:
		
		1) YCbCr 1x1 1x1 1x1
		2) YCbCr 2x1 1x1 1x1
		3) YCbCr 2x2 1x1 1x1
		
  * Huffman encoding only.
  
  * No saving/encoding supported, this library is for loading/decoding only.
    
Error handling
--------------
  C++ exceptions are used to indicate decoding failures. Only a single 
  exception class can be thrown: ttjpg::DecodeError. 
  
  To gain more information about the failure use the getError() method.
  
  For a list of possible error codes consult TinyTinyJPG.h.
  
Usage
-----
  This library is very easy to integrate in any C++ application. The task
  can be summarized as follows:
  
  1) Create an instance of ttjpg::Decoder which is the main JPEG decoding 
     object. If you're going to decode multiple JPEGs it's useful to use the 
     same decoding object, since it does some initialization and 
     precalculations.
     
     Note that the entire library resides in the ttjpg namespace.
  
  2) You need to tell ttjpg where it's supposed to read the binary JPEG data.
     This is done by implementing the ttjpg::InputStream class. You'll need to
     implement the following virtual methods:
     
       virtual int readByte(void);       
       virtual bool isEndOfStream(void);
       virtual void skip(int nBytes);

     readByte() should read a single byte and return the value.       
     isEndStream() should return true if a previous call to readByte() 
     attempted to read past the end of the stream (file), false otherwise.
     skip() should skip the next nBytes bytes of the stream.
     
     See the StdioInput class defined in test/Test.cpp for an example.
     
     Note that you should probably do some memory buffering, mapping readByte()
     directly to, for example, fgetc() can be very slow due to multithreading
     locks.
  
  3) When you got your very own ttjpg::InputStream-implementing stream
     up and running, it's just a matter of passing the stream to the 
     JPEG decoder object. The decoder object has two methods:
     
       void readImage(InputStream *pIn,RGBImage *pRGB);    		
       
         This method read a JPEG image from the given stream and stores it in
         the supplied ttjpg::RGBImage struct. If a previous RGB image is still
         stored in it, it will be released first.
         
         If the method is succesful it returns and you'll find the following in
         the RGBImage struct:
         
           nWidth      The width of the image in pixels.
           nHeight     The height of the image in pixels.
           pc          A pointer to the raw RGB image (24-bit interleaved 
                       format, 8 bits per component). Note that this will be 
                       released automatically when the ttjpg::RGBImage runs out
                       of scope.
         
         If the decoding fails a DecodeError exception is thrown.
       
       bool readImageInfo(InputStream *pIn,RGBImage *pRGB);
		 
		 Acts much like readImage(), except it will stop before actually 
		 decoding the image. The width and height members of the supplied
		 ttjpg:RGBImage struct will be set, pc will be NULL.

Build options and dependencies
------------------------------
  No special considerations needed, just add TinyTinyJPG.cpp to your project/
  makefile/whatever. 
  
  To enable SSE2 optimizations (Visual C++ only), #define EMMINTRIN. 
  
  Visual Studio 2008 solution is supplied in the vs2008 directory which will
  build the unit++ based test suite.
  
License
-------
  The source code is available under the BSD license (see LICENSE.txt for 
  details). 
