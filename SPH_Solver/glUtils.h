#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "GL/glew.h"
#include "GL/freeglut.h"

#include "geometry.h"

// Creating Shaders


void ExitOnGLError(const char* error_message);
GLuint LoadShader(const char* filename, GLenum shader_type);