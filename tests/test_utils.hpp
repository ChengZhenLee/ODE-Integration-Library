#pragma once


#include <iomanip>
#include <iostream>


namespace TestUtils {

inline void printHeader(const std::string& title) {
    std::cout << "\n-- " << title << " " << std::string(std::max<int>(0, 40 - static_cast<int>(title.size())), '-') << "\n";
}

template <typename T>
inline void printRow(const std::string& label, const T& value) {
    std::cout << "  " << std::left << std::setw(22) << label
               << std::right << std::scientific << std::setprecision(6) << value << "\n";
}

inline void printPassed(const std::string& name) {
    std::cout << "[PASSED] " << name << "\n";
};

}
