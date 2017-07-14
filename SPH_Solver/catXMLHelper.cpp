#include "CatToolBox.h"

using namespace tinyxml2;


cfloat3 QueryFloat3(XMLElement* element) {
	cfloat3 tmp(0,0,0);
	if(element != NULL){
		const char* text = element->GetText();
		sscanf(text, "%f,%f,%f", &tmp.x, &tmp.y, &tmp.z);
	}
	return tmp;
}