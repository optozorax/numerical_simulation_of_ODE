#pragma once

#include <vector>
#include <iostream>

typedef std::vector<double> vec;

inline vec operator+(const vec& a, const vec& b) {
	#ifdef _DEBUG
	if (a.size() != b.size())
		throw std::exception();
	#endif
	vec result = a;
	for (int i = 0; i < b.size(); i++)
		result[i] += b[i];
	return result;
}

inline vec operator-(const vec& a, const vec& b) {
	#ifdef _DEBUG
	if (a.size() != b.size())
		throw std::exception();
	#endif
	vec result = a;
	for (int i = 0; i < b.size(); i++)
		result[i] -= b[i];
	return result;
}

inline vec operator*(const vec& a, double b) {
	vec result = a;
	for (int i = 0; i < result.size(); i++)
		result[i] *= b;
	return result;
}

inline vec operator*(double b, const vec& a) {
	return operator*(a, b);
}

inline std::ostream& operator<<(std::ostream& out, const vec& v) {
	for (int i = 0; i < v.size()-1; ++i)
		out << v[i] << "\t";
	out << v.back();
	return out;
}