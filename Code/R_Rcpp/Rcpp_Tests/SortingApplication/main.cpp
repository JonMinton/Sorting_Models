/* 
 * File:   main.cpp
 * Author: Jon Minton
 *
 * Created on 20 January 2014, 22:35
 */


#include <istream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdio>
#include <iostream>


using namespace std;


int main() {
    
    ifstream file_medval, file_perc_renter, file_pop_data;
    vector <int> medval;
    vector <long> perc_renter;
    vector <int> pop_data;
    int number, i;
    long numlong;
    
    
    file_medval.open("X:/RcppTest/SortingApplication/inputData/medval.txt", ios::in);
    
    if(file_medval.is_open())
    {
        while (file_medval>>number)
        {
            medval.push_back(number);
            file_medval.get();
        }
    }
    file_medval.close();

    file_perc_renter.open("X:/RcppTest/SortingApplication/inputData/perc_renter.txt", ios::in);
    
    if(file_perc_renter.is_open())
    {
        while (file_perc_renter>>numlong)
        {
            perc_renter.push_back(numlong);
            file_perc_renter.get();
        }
    }
    file_perc_renter.close();

    
    file_pop_data.open("X:/RcppTest/SortingApplication/inputData/pop_data.txt", ios::in);
    
    if(file_pop_data.is_open())
    {
        while (file_pop_data>>number)
        {
            pop_data.push_back(number);
            file_pop_data.get();
        }
    }
    file_pop_data.close();
    
    cout << medval.size() << "\t" << perc_renter.size() << "\t" << pop_data.size() << "\n";    
    
    return 0;
}

