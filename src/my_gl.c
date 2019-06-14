#include "my_gl.h"
#include <stdio.h>
#if 0
void __stdcall my_glVertex3fv (const GLfloat *v)
{
	glVertex3fv( v );
	//printf("glVertex3fv %f %f %f\n", v[0], v[1], v[2]);
}
void __stdcall my_glVertex3dv (const GLdouble *v)
{
	glVertex3dv( v );
	//printf("glVertex3dv %f %f %f\n", v[0], v[1], v[2]);
}

void __stdcall my_glNormal3fv (const GLfloat *v)
{
	glNormal3fv( v );
	//printf("glNormal3fv %f %f %f\n", v[0], v[1], v[2]);
}

void __stdcall my_glNormal3dv (const GLdouble *v)
{
	glNormal3dv( v );
	//printf("glNormal3dv %f %f %f\n", v[0], v[1], v[2]);
}

void __stdcall my_glVertex2fv (const GLfloat *v)
{
	glVertex2fv( v );
	//printf("glVertex2fv %f %f\n", v[0], v[1]);
}

void __stdcall my_glVertex2dv (const GLdouble *v)
{
	glVertex2dv( v );
}

typedef void (WINAPI* PFNGLMULTITEXCOORD2FVPROC) (GLenum target, const GLfloat *v);
typedef void (WINAPI* PFNGLMULTITEXCOORD2DVPROC) (GLenum target, const GLdouble *v);

static PFNGLMULTITEXCOORD2FVPROC _my_glMultiTexCoord2fv;
static PFNGLMULTITEXCOORD2DVPROC _my_glMultiTexCoord2dv;

static unsigned int gl_version = 0;

void my_GetGLVersion()
{
	if ( gl_version == 0 )
	{
		/* check version and extensions */
		int uiMajor, uiMinor;
		unsigned char v12,v13,v14,v15,v20,v21,v30;
	
		sscanf((char*)glGetString(GL_VERSION), "%d.%d", &uiMajor, &uiMinor);
	
		gl_version = (uiMajor<<8) | uiMinor;
	}
}

void __stdcall my_glMultiTexCoord2fv( GLenum target, const GLfloat *v )
{
	unsigned char v12,v13,v14,v15,v20,v21,v30;
	
	my_GetGLVersion();

	v12 = gl_version >= 0x102;
	v13 = gl_version >= 0x103;
	v14 = gl_version >= 0x104;
	v15 = gl_version >= 0x105;
	v20 = gl_version >= 0x200;
	v21 = gl_version >= 0x201;
	v30 = gl_version >= 0x300;

	_my_glMultiTexCoord2fv = (PFNGLMULTITEXCOORD2FVPROC)
			wglGetProcAddress(v13 ? "glMultiTexCoord2fv" : "glMultiTexCoord2fvARB");

	_my_glMultiTexCoord2fv(target, v);
	//printf("glMultiTexCoord2fv %f %f\n", v[0], v[1]);
}


void __stdcall my_glMultiTexCoord2dv( GLenum target, const GLdouble *v )
{
	unsigned char v12,v13,v14,v15,v20,v21,v30;
	
	my_GetGLVersion();

	v12 = gl_version >= 0x102;
	v13 = gl_version >= 0x103;
	v14 = gl_version >= 0x104;
	v15 = gl_version >= 0x105;
	v20 = gl_version >= 0x200;
	v21 = gl_version >= 0x201;
	v30 = gl_version >= 0x300;

	_my_glMultiTexCoord2dv = (PFNGLMULTITEXCOORD2DVPROC)
			wglGetProcAddress(v13 ? "glMultiTexCoord2dv" : "glMultiTexCoord2dvARB");

	_my_glMultiTexCoord2dv(target, v);
	//printf("glMultiTexCoord2dv %f %f\n", v[0], v[1]);
}
#endif