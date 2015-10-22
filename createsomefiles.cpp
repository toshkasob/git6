/* #include <vector> #include <iostream> #include <fstream> using namespace
   std; int main() { fstream file; file.open(

   } */
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
int main()
{
	char q[255] = "test";
	char buffer[33];
	for (int i = 0; i < 10; i++)
	{
		q[4] = 0;
		itoa(i, buffer, 10);
		strcat(q, buffer);
		strcat(q, ".txt");
		ofstream out(q);
		out << "test";
	}
	system("pause");
	return 0;
}
// http://www.cyberforum.ru/cpp-beginners/thread254122.html