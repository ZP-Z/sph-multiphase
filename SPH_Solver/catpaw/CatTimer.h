#pragma once

#include<Windows.h>
#include <time.h>
#include <stdio.h>

struct cTime{
	SYSTEMTIME st;
	
	cTime(){
		GetSystemTime(&st);
	}

	inline cTime& update(){
		GetSystemTime(&st);
		return *this;
	}

	inline int dailyMilliSec(){
		return st.wHour*3600000+st.wMinute*60000+st.wSecond*1000 + st.wMilliseconds;
	}
	inline int operator- (cTime& b){
		return dailyMilliSec() - b.dailyMilliSec();
	}
	inline void printTime(){
		printf("%hu:%hu:%hu:%hu",st.wHour,st.wMinute,st.wSecond,st.wMilliseconds);
	}
};