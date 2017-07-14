
#pragma once

/*#define TEX_SIZE		2048	
#define LIGHT_NEAR		0.5
#define LIGHT_FAR		300.0
#define DEGtoRAD		(3.141592/180.0)*/

typedef unsigned int			uint;
typedef unsigned short int		ushort;

/*#define CLRVAL			uint
#define COLOR(r,g,b)	( (uint(r*255.0f)<<24) | (uint(g*255.0f)<<16) | (uint(b*255.0f)<<8) )
#define COLORA(r,g,b,a)	( (uint(a*255.0f)<<24) | (uint(b*255.0f)<<16) | (uint(g*255.0f)<<8) | uint(r*255.0f) )
#define ALPH(c)			(float((c>>24) & 0xFF)/255.0)
#define BLUE(c)			(float((c>>16) & 0xFF)/255.0)
#define GRN(c)			(float((c>>8)  & 0xFF)/255.0)
#define RED(c)			(float( c      & 0xFF)/255.0)*/

#ifndef EMIT_BUF_RATIO
#define EMIT_BUF_RATIO 2
#endif

#define MAX_FLUIDNUM 3