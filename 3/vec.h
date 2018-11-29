#pragma once

#include <vector>

typedef std::vector<double> vec;

vec operator+(const vec& a, const vec& b) {
	#ifdef _DEBUG
	if (a.size() != b.size())
		throw std::exception();
	#endif
	vec result = a;
	for (int i = 0; i < b.size(); i++)
		result[i] += b[i];
	return result;
}

vec operator-(const vec& a, const vec& b) {
	#ifdef _DEBUG
	if (a.size() != b.size())
		throw std::exception();
	#endif
	vec result = a;
	for (int i = 0; i < b.size(); i++)
		result[i] -= b[i];
	return result;
}

vec operator*(const vec& a, double b) {
	vec result = a;
	for (int i = 0; i < result.size(); i++)
		result[i] *= b;
	return result;
}

vec operator*(double b, const vec& a) {
	return operator*(a, b);
}