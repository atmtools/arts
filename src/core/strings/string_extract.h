#pragma once

#include <string>

void extract(float& x, std::string& line, std::size_t n);
void extract(double& x, std::string& line, std::size_t n);
void extract(char& x, std::string& line, std::size_t n);
void extract(int& x, std::string& line, std::size_t n);
void extract(long int& x, std::string& line, std::size_t n);
void extract(long long int& x, std::string& line, std::size_t n);
void extract(unsigned char& x, std::string& line, std::size_t n);
void extract(unsigned int& x, std::string& line, std::size_t n);
void extract(unsigned long int& x, std::string& line, std::size_t n);
void extract(unsigned long long int& x, std::string& line, std::size_t n);
