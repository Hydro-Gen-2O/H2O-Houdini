#ifndef COMMON_DEF
	#define COMMON_DEF

	// Global defs
	#define DEGtoRAD		(3.141592/180.0)

	#include <windows.h>
	typedef unsigned int		uint;

	#define COLOR(r,g,b)	( (DWORD(r*255.0f)<<24) | (DWORD(g*255.0f)<<16) | (DWORD(b*255.0f)<<8) )
	#define COLORA(r,g,b,a)	( (DWORD(r*255.0f)<<24) | (DWORD(g*255.0f)<<16) | (DWORD(b*255.0f)<<8) | DWORD(a*255.0f) )
	#define RED(c)			(float((c>>24) & 0xFF)/255.0)
	#define GRN(c)			(float((c>>16) & 0xFF)/255.0)
	#define BLUE(c)			(float((c>>8) & 0xFF)/255.0)
	#define ALPH(c)			(float(c & 0xFF)/255.0)
#endif