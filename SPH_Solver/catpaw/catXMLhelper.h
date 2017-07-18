#pragma once
#include "tinyxml2.h"
#include <string>


using namespace tinyxml2;

cfloat3 QueryFloat3(XMLElement* element);
float XMLGetFloat(const char* name);
int XMLGetInt(const char* name);
cfloat3& XMLGetFloat3(const char* name);
void XMLGetFloatN(float* buffer, int size, const char* name);