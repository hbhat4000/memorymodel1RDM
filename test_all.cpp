#include <Eigen/Dense>
#include <iostream>

int main() {
    auto a = Eigen::placeholders::all;
    std::cout << "Eigen::all successfully compiled!" << std::endl;
    return 0;
}
