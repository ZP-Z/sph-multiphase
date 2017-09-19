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


void Tinyxml_Reader::Use(XMLElement* base) {
	basenode = base;
	if (!base) {
		printf("WARNING: NULL Node Used.\n");
	}
}

const char* Tinyxml_Reader::GetText(const char* name) {
	float tmp=0;
	if (basenode!=NULL) {
		XMLElement* attr = basenode->FirstChildElement(name);
		if (attr)
			return attr->GetText();
	}
	return NULL;
}

float Tinyxml_Reader::GetFloat(const char* name) {
	float tmp=0;
	if (basenode!=NULL) {
		XMLElement* attr = basenode->FirstChildElement(name);
		if (attr)
			attr->QueryFloatText(&tmp);
	}
	return tmp;
}

int Tinyxml_Reader::GetInt(const char* name) {
	int tmp=0;
	if (basenode!=NULL) {
		basenode->FirstChildElement(name)->QueryIntText(&tmp);
	}
	return tmp;
}

cfloat3& Tinyxml_Reader::GetFloat3(const char* name) {
	if (basenode!=NULL) {
		return QueryFloat3(basenode->FirstChildElement(name));
	}
	else
		return cfloat3(0, 0, 0);
}



void Tinyxml_Reader::GetFloatN(float* buffer, int size, const char* name) {
	if (basenode==NULL)
		return;
	const char* str = basenode->FirstChildElement(name)->GetText();
	std::string fmstr = "";
	for (int i=0; i<size-1; i++)
		fmstr = fmstr+"%f,";
	fmstr = fmstr+"%f";
	sscanf(str, fmstr.c_str(), &buffer[0], &buffer[1], &buffer[2]);
}