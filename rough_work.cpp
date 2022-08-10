#include <iostream>
#include <fstream>
using namespace std;

int main()
{
    fstream myFile;

    myFile.open("Text.txt", ios::out);
    if (myFile.is_open())
    {
        myFile << "Hello!\n";
        myFile << "This is 2nd line.\n";
        myFile.close();
    }
    

}