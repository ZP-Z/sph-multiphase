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

XMLElement* pele=NULL;
float XMLGetFloat(const char* name) {
	float tmp=0;
	if (pele!=NULL) {
		XMLElement* attr = pele->FirstChildElement(name);
		if (attr)
			attr->QueryFloatText(&tmp);
	}
	return tmp;
}

int XMLGetInt(const char* name) {
	int tmp=0;
	if (pele!=NULL) {
		pele->FirstChildElement(name)->QueryIntText(&tmp);
	}
	return tmp;
}

cfloat3& XMLGetFloat3(const char* name) {
	if (pele!=NULL) {
		return QueryFloat3(pele->FirstChildElement(name));
	}
	else
		return cfloat3(0, 0, 0);
}



void XMLGetFloatN(float* buffer, int size, const char* name) {
	if (pele==NULL)
		return;
	const char* str = pele->FirstChildElement(name)->GetText();
	std::string fmstr = "";
	for (int i=0; i<size-1; i++)
		fmstr = fmstr+"%f,";
	fmstr = fmstr+"%f";
	sscanf(str, fmstr.c_str(), &buffer[0], &buffer[1], &buffer[2]);
}
