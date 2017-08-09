#include<iostream>
#include<fstream>
#include<math.h>
#include<stdio.h>
#include<cstdlib>
#include<cmath>
#include<string>
#include<sstream>

#define pi 3.14159265
#define NUMSIMULATIONS 5000000

using namespace std;



int main()
{
	
	string FirstName = "m2m46_red_";
	string LastName = "_";
	for(int age = 0; age < 10; age++)
	{
		ostringstream convert;
		ostringstream convert2;
		convert << age;
		convert2 << NUMSIMULATIONS;
		string FullName = FirstName + convert.str() + LastName + convert2.str();
		ofstream myfile3((FullName+".txt").c_str());
		if(myfile3.is_open())
		{
			myfile3 << FullName;
		}
	}
	
	
	
	
	return 0;
}