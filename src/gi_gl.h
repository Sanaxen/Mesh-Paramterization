/*
 *  OpenGI: Library for Parameterization and Geometry Image creation
 *  Copyright (C) 2008-2011  Christian Rau
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published 
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library. If not, see <http://www.gnu.org/licenses/>.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Christian Rau
 *
 *     rauy@users.sourceforge.net
 */

/** \internal
 *  \file
 *  \brief Declaration of structures and functions for communication with OpenGL.
 */

#ifndef __GI_GL_H__
#define __GI_GL_H__

#if HAVE_CONFIG_H
	#include <config.h>
#endif
#include <GI/gi.h>

#include "gi_container.h"
#include "gi_mesh.h"
#include "gi_cutter.h"

//#ifdef _WIN32
//	#include <windows.h>
//	#define GI_GL_PROC_ADDRESS(str)		wglGetProcAddress(str)
//#else
//	#include <GL/glx.h>
//	#define GI_GL_PROC_ADDRESS(str)		glXGetProcAddress(str)
//#endif
//
//#include <GL/gl.h>
//#include <GL/glext.h>

#include "my_types.h"

/*************************************************************************/
/* Macros */

#define GI_GL_SEMANTIC_BASE				GI_GL_VERTEX
#define GI_GL_SEMANTIC_END				GI_GL_VERTEX_ATTRIB
#define GI_GL_SEMANTIC_COUNT			(GI_GL_SEMANTIC_END-GI_GL_SEMANTIC_BASE+1)

#define GI_NORMALIZATION_CUBE_SIZE		64

#define GI_SAMPLE_SHADER_DEFAULT		0
#define GI_SAMPLE_SHADER_TEXTURE1D		1
#define GI_SAMPLE_SHADER_TEXTURE2D		2
#define GI_SAMPLE_SHADER_TEXTURE3D		3
#define GI_SAMPLE_SHADER_TEXTURECUBE	4
#define GI_SAMPLE_SHADER_NORMALIZED		5
#define GI_SAMPLE_SHADER_PACKED			6
#define GI_NUM_SAMPLERS					7


/*************************************************************************/
/* Typedefs */

/** \internal
 *  \brief Function for setting single-precision GL attributes.
 *  \ingroup opengl
 */
//typedef void (APIENTRY *GIGLattribffunc)(const GLfloat*);

/** \internal
 *  \brief Function for setting double-precision GL attributes.
 *  \ingroup opengl
 */
//typedef void (APIENTRY *GIGLattribdfunc)(const GLdouble*);

/** \internal
 *  \brief Function for setting single-precision 
 *  multitexture coordinates and vertex attributes.
 *  \ingroup opengl
 */
//typedef void (APIENTRY *GIGLidxattribffunc)(GLuint, const GLfloat*);

/** \internal
 *  \brief Function for setting double-precision 
 *  multitexture coordinates and vertex attributes.
 *  \ingroup opengl
 */
//typedef void (APIENTRY *GIGLidxattribdfunc)(GLuint, const GLdouble*);


/*************************************************************************/
/* Structures */

/** \internal
 *  \brief Render state.
 *  \ingroup opengl
 */
typedef struct _GIRenderer
{
	struct _GIContext	*context;							/**< Context this renderer belongs to */
	GIenum				attrib_semantic[GI_ATTRIB_COUNT];	/**< Render semantics of attributes. */
	GIint				attrib_channel[GI_ATTRIB_COUNT];	/**< Render channels of attributes. */
	GIsizei				render_res[2];						/**< Resolution for GIM rendering. */
	GIsizei				gim_cache_size;						/**< Size of geometry image cache. */
	GIenum				image_domain[GI_ATTRIB_COUNT];		/**< TexCoord domains for images. */
	GIbitfield			gim_flags;							/**< GIM render flags. */
} GIRenderer;

/** \internal
 *  \brief Attribute render information.
 *  \details This structure represents the information about an attribute to render.
 *  \ingroup opengl
 */
typedef struct _GIRenderAttrib
{
	GIuint	attrib;								/**< Attribute channel. */
	GIint	index;								/**< Render index (attrib location or texture unit). */
	GIvoid	*func;								/**< GL attribute function. */
} GIRenderAttrib;

/** \internal
 *  \brief GLSL shader.
 *  \details This structure represents a shader program in OpenGL shading language.
 *  \ingroup opengl
 */
typedef struct _GIGLShader
{
	GIuint	program;							/**< OpenGL program object. */
	GIuint	v_shader;							/**< Shader object for vertex shader. */
	GIuint	g_shader;							/**< Shader object for geometry shader. */
	GIuint	f_shader;							/**< Shader object for fragment shader. */
	GIHash	attrib_locs;						/**< Table of attribute locations. */
	GIHash	uniform_locs;						/**< Table of uniform locations. */
} GIGLShader;

/** \internal
 *  \brief Geometry image cache element.
 *  \ingroup opengl
 */
typedef struct _GIGLGIMCache
{
	GIuint					width_height;		/**< Packed image width and height. */
	GIuint					res_uv;				/**< Packed rendering resolution. */
	GIuint					buffer;				/**< GL buffer object. */
	struct _GIGLGIMCache	*next;				/**< Next cache element in list. */
	struct _GIGLGIMCache	*prev;				/**< Previous cache element in list. */
} GIGLGIMCache;

/** \internal
 *  \brief OpenGL extension manager.
 *  \details This structure manages OpenGL extensions by supplying boolean flags and function pointers.
 *  \ingroup opengl
 */
typedef struct _GIGLManager
{
	GIuint								gl_version;						/**< OpenGL version (LSB minor, next byte major). */
	GIuint								glsl_version;					/**< GLSL version (LSB minor, next byte major). */
	GIboolean							texture_3d;						/**< 3D textures supported. */
	GIboolean							multitexture;					/**< Multitexturing supported. */
	GIboolean							texture_cube_map;				/**< Cube mapping supported. */
	GIboolean							non_power_of_2;					/**< Non power of two textures supported. */
	GIboolean							vbo;							/**< Vertex buffer objects supported. */
	GIboolean							pbo;							/**< Pixel buffer objects supported. */
	GIboolean							glsl;							/**< GLSL shaders supported. */
	GIboolean							vertex_texture;					/**< Vertex textureing supported. */
	GIboolean							geometry_shader;				/**< Geometry shader supported. */
	GIboolean							fbo;							/**< Framebuffer objects supported. */
	GIboolean							color_buffer_float;				/**< Floating point color buffer supported. */
	GIboolean							half_float_pixel;				/**< Floating point pixel data supported. */
	GIboolean							half_float_vertex;				/**< Floating point vertex data supported. */
	GIboolean							texture_float;					/**< Floating point textures supported. */
	GIboolean							texture_rg;						/**< RG textures supporte. */
	GIboolean							draw_buffers;					/**< Multiple draw buffers supported. */
	GIint								max_vertices;					/**< Maximal number of vertices for range elements. */
	GIint								max_indices;					/**< Maximal number of indices for range elements. */
	GIint								max_texunits;					/**< Maximal number of texture units. */
	GIint								max_attachments;				/**< Maximum number of color attachments for FBOs. */
	GIint								max_draw_buffers;				/**< Maximum number of draw buffers. */
	GIuint								norm_cubemap[3];				/**< Normalization cube maps. */
	GIuint								gim_buffer;						/**< Preallocated buffer for geometry data. */
	GIuint								gim_cache_size;					/**< Size of geometry image cache. */
	GIGLGIMCache						*texCoord_cache;				/**< Cache for texCoord data. */
	GIGLGIMCache						*index_cache;					/**< Cache for index data. */
	//GIGLShader							*gim_sampler[GI_NUM_SAMPLERS];	/**< Shaders for geometry image sampling. */
	//GIGLShader							*gim_renderer[4];				/**< Shaders for geometry image rendering. */
	//PFNGLDRAWRANGEELEMENTSPROC			_glDrawRangeElements;			/**< Pointer to glDrawRangeElements function. */
	//PFNGLACTIVETEXTUREPROC				_glActiveTexture;				/**< Pointer to glActiveTexture function. */
	//PFNGLCLIENTACTIVETEXTUREPROC		_glClientActiveTexture;			/**< Pointer to glClientActiveTexture function. */
	//PFNGLMULTITEXCOORD1IPROC			_glMultiTexCoord1i;				/**< Pointer to glMultiTexCoord1i function. */
	//PFNGLMULTITEXCOORD1FPROC			_glMultiTexCoord1f;				/**< Pointer to glMultiTexCoord1f function. */
	//PFNGLMULTITEXCOORD1FVPROC			_glMultiTexCoord1fv;			/**< Pointer to glMultiTexCoord1fv function. */
	//PFNGLMULTITEXCOORD2FVPROC			_glMultiTexCoord2fv;			/**< Pointer to glMultiTexCoord2fv function. */
	//PFNGLMULTITEXCOORD3FVPROC			_glMultiTexCoord3fv;			/**< Pointer to glMultiTexCoord3fv function. */
	//PFNGLMULTITEXCOORD4FVPROC			_glMultiTexCoord4fv;			/**< Pointer to glMultiTexCoord4fv function. */
	//PFNGLMULTITEXCOORD1DPROC			_glMultiTexCoord1d;				/**< Pointer to glMultiTexCoord1d function. */
	//PFNGLMULTITEXCOORD1DVPROC			_glMultiTexCoord1dv;			/**< Pointer to glMultiTexCoord1fv function. */
	//PFNGLMULTITEXCOORD2DVPROC			_glMultiTexCoord2dv;			/**< Pointer to glMultiTexCoord2fv function. */
	//PFNGLMULTITEXCOORD3DVPROC			_glMultiTexCoord3dv;			/**< Pointer to glMultiTexCoord3fv function. */
	//PFNGLMULTITEXCOORD4DVPROC			_glMultiTexCoord4dv;			/**< Pointer to glMultiTexCoord4fv function. */
	//PFNGLFOGCOORDFVPROC					_glFogCoordfv;					/**< Pointer to glFogCoordfv function. */
	//PFNGLFOGCOORDDVPROC					_glFogCoorddv;					/**< Pointer to glFogCoorddv function. */
	//PFNGLSECONDARYCOLOR3FVPROC			_glSecondaryColor3fv;			/**< Pointer to glSecondaryColor3fv function. */
	//PFNGLSECONDARYCOLOR3DVPROC			_glSecondaryColor3dv;			/**< Pointer to glSecondaryColor3dv function. */
	//PFNGLBINDBUFFERPROC					_glBindBuffer;					/**< Pointer to glBindBuffer function. */
	//PFNGLBUFFERDATAPROC					_glBufferData;					/**< Pointer to glBufferData function. */
	//PFNGLBUFFERSUBDATAPROC				_glBufferSubData;				/**< Pointer to glBufferSubData function. */
	//PFNGLDELETEBUFFERSPROC				_glDeleteBuffers;				/**< Pointer to glDeleteBuffers function. */
	//PFNGLGENBUFFERSPROC					_glGenBuffers;					/**< Pointer to glGenBuffers function. */
	//PFNGLGETBUFFERSUBDATAPROC			_glGetBufferSubData;			/**< Pointer to glGetBufferSubData function. */
	//PFNGLISBUFFERPROC					_glIsBuffer;					/**< Pointer to glIsBuffer function. */
	//PFNGLMAPBUFFERPROC					_glMapBuffer;					/**< Pointer to glMapBuffer function. */
	//PFNGLUNMAPBUFFERPROC				_glUnmapBuffer;					/**< Pointer to glUnmapBuffer function. */
	//PFNGLATTACHSHADERPROC				_glAttachShader;				/**< Pointer to glAttachShader function. */
	//PFNGLCOMPILESHADERPROC				_glCompileShader;				/**< Pointer to glCompileShader function. */
	//PFNGLCREATEPROGRAMPROC				_glCreateProgram;				/**< Pointer to glCreateProgram function. */
	//PFNGLCREATESHADERPROC				_glCreateShader;				/**< Pointer to glCreateShader function. */
	//PFNGLDELETEPROGRAMPROC				_glDeleteProgram;				/**< Pointer to glDeleteProgram function. */
	//PFNGLDELETESHADERPROC				_glDeleteShader;				/**< Pointer to glDeleteShader function. */
	//PFNGLDETACHSHADERPROC				_glDetachShader;				/**< Pointer to glDetachShader function. */
	//PFNGLDISABLEVERTEXATTRIBARRAYPROC	_glDisableVertexAttribArray;	/**< Pointer to glDisableVertexAttribArray function. */
	//PFNGLENABLEVERTEXATTRIBARRAYPROC	_glEnableVertexAttribArray;		/**< Pointer to glEnableVertexAttribArray function. */
	//PFNGLGETACTIVEATTRIBPROC			_glGetActiveAttrib;				/**< Pointer to glGetActiveAttrib function. */
	//PFNGLGETACTIVEUNIFORMPROC			_glGetActiveUniform;			/**< Pointer to glGetActiveUniform function. */
	//PFNGLGETATTRIBLOCATIONPROC			_glGetAttribLocation;			/**< Pointer to glGetAttribLocation function. */
	//PFNGLGETPROGRAMINFOLOGPROC			_glGetProgramInfoLog;			/**< Pointer to glGetProgramInfoLog function. */
	//PFNGLGETPROGRAMIVPROC				_glGetProgramiv;				/**< Pointer to glGetProgramiv function. */
	//PFNGLGETSHADERINFOLOGPROC			_glGetShaderInfoLog;			/**< Pointer to glGetShaderInfoLog function. */
	//PFNGLGETSHADERIVPROC				_glGetShaderiv;					/**< Pointer to glGetShaderiv function. */
	//PFNGLGETUNIFORMLOCATIONPROC			_glGetUniformLocation;			/**< Pointer to glGetUniformLocation function. */
	//PFNGLLINKPROGRAMPROC				_glLinkProgram;					/**< Pointer to glLinkProgram function. */
	//PFNGLSHADERSOURCEPROC				_glShaderSource;				/**< Pointer to glShaderSource function. */
	//PFNGLUNIFORM1IPROC					_glUniform1i;					/**< Pointer to glUniform1i function. */
	//PFNGLUNIFORM2IPROC					_glUniform2i;					/**< Pointer to glUniform2i function. */
	//PFNGLUNIFORM3IPROC					_glUniform3i;					/**< Pointer to glUniform3i function. */
	//PFNGLUNIFORM4IPROC					_glUniform4i;					/**< Pointer to glUniform4i function. */
	//PFNGLUNIFORM1IVPROC					_glUniform1iv;					/**< Pointer to glUniform4iv function. */
	//PFNGLUNIFORM4IVPROC					_glUniform4iv;					/**< Pointer to glUniform4iv function. */
	//PFNGLUNIFORM2FPROC					_glUniform2f;					/**< Pointer to glUniform2f function. */
	//PFNGLUNIFORM2FVPROC					_glUniform2fv;					/**< Pointer to glUniform2fv function. */
	//PFNGLUNIFORM4FVPROC					_glUniform4fv;					/**< Pointer to glUniform4fv function. */
	//PFNGLUSEPROGRAMPROC					_glUseProgram;					/**< Pointer to glUseProgram function. */
	//PFNGLVALIDATEPROGRAMPROC			_glValidateProgram;				/**< Pointer to glValidateProgram function. */
	//PFNGLVERTEXATTRIB1FVPROC			_glVertexAttrib1fv;				/**< Pointer to glVertexAttrib1fv function. */
	//PFNGLVERTEXATTRIB2FVPROC			_glVertexAttrib2fv;				/**< Pointer to glVertexAttrib2fv function. */
	//PFNGLVERTEXATTRIB3FVPROC			_glVertexAttrib3fv;				/**< Pointer to glVertexAttrib3fv function. */
	//PFNGLVERTEXATTRIB4FVPROC			_glVertexAttrib4fv;				/**< Pointer to glVertexAttrib4fv function. */
	//PFNGLVERTEXATTRIB1DVPROC			_glVertexAttrib1dv;				/**< Pointer to glVertexAttrib1dv function. */
	//PFNGLVERTEXATTRIB2DVPROC			_glVertexAttrib2dv;				/**< Pointer to glVertexAttrib2dv function. */
	//PFNGLVERTEXATTRIB3DVPROC			_glVertexAttrib3dv;				/**< Pointer to glVertexAttrib3dv function. */
	//PFNGLVERTEXATTRIB4DVPROC			_glVertexAttrib4dv;				/**< Pointer to glVertexAttrib4dv function. */
	//PFNGLVERTEXATTRIBPOINTERPROC		_glVertexAttribPointer;			/**< Pointer to glVertexAttribPointer function. */
	//PFNGLGETHANDLEARBPROC				_glGetHandle;					/**< Pointer to glGetHandle function. */
	//PFNGLPROGRAMPARAMETERIEXTPROC		_glProgramParameteri;			/**< Pointer to glProgramParameteri function. */
	//PFNGLBINDFRAMEBUFFEREXTPROC			_glBindFramebuffer;				/**< Pointer to glBindFramebuffer function. */
	//PFNGLBINDRENDERBUFFEREXTPROC		_glBindRenderbuffer;			/**< Pointer to glBindRenderbuffer function. */
	//PFNGLCHECKFRAMEBUFFERSTATUSEXTPROC	_glCheckFramebufferStatus;		/**< Pointer to glCheckFramebufferStatus function. */
	//PFNGLDELETEFRAMEBUFFERSEXTPROC		_glDeleteFramebuffers;			/**< Pointer to glFramebufferRenderbuffer function. */
	//PFNGLDELETERENDERBUFFERSEXTPROC		_glDeleteRenderbuffers;			/**< Pointer to glDeleteRenderbuffers function. */
	//PFNGLFRAMEBUFFERRENDERBUFFEREXTPROC	_glFramebufferRenderbuffer;		/**< Pointer to glFramebufferRenderbuffer function. */
	//PFNGLFRAMEBUFFERTEXTURE2DEXTPROC	_glFramebufferTexture2D;		/**< Pointer to glFramebufferTexture2D function. */
	//PFNGLGENFRAMEBUFFERSEXTPROC			_glGenFramebuffers;				/**< Pointer to glGenFramebuffers function. */
	//PFNGLGENRENDERBUFFERSEXTPROC		_glGenRenderbuffers;			/**< Pointer to glGenRenderbuffers function. */
	//PFNGLRENDERBUFFERSTORAGEEXTPROC		_glRenderbufferStorage;			/**< Pointer to glRenderbufferStorage function. */
	//PFNGLDRAWBUFFERSPROC				_glDrawBuffers;					/**< Pointer to glDrawBuffers function. */
	//PFNGLCLAMPCOLORPROC					_glClampColor;					/**< Pointer to glClampColor function. */
} GIGLManager;


/*************************************************************************/
/* Functions */

/** \name Renderer methods
 *  \{
 */
void GIRenderer_construct(GIRenderer *renderer, struct _GIContext *context);
GIboolean GIRenderer_setup_mesh_rendering(GIRenderer *renderer, 
	GIRenderAttrib *attribs, GIuint *num_attribs, GIPatch **pstart, GIPatch **pend);
/** \} */

/** \name Shader methods
 *  \{
 */
void GIGLShader_construct(GIGLShader *shader, GIGLManager *gl, 
	const GIchar *vsource, const GIchar *gsource, const GIchar *fsource, 
	GIenum g_in, GIenum g_out, GIuint num_g_out, GIuint num_defines, ...);
void GIGLShader_destruct(GIGLShader *shader, GIGLManager *gl);
/** \} */

/** name GIM cache methods
 *  \{
 */
void GIGLGIMCache_destruct(GIGLGIMCache *cache, GIGLManager *gl);
/** \} */

/** \name GL manager methods
 *  \{
 */
void GIGLManager_init(GIGLManager *mgr);
void GIGLManager_destruct(GIGLManager *mgr);
GIuint GIGLManager_normalization_cubemap(GIGLManager *mgr, GIenum type);
/** \} */


#endif