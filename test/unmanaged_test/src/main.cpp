#include <stdlib.h>
#include <iostream>
#include "__float128.h"

// Utility
__float128 fRand(__float128 min, __float128 max) {
    __float128 f = (__float128)rand() / RAND_MAX;
    return min + f * (max-min);
}
int iRand(int min, int max) {
    double f = (double)rand() / RAND_MAX;
    return min + (int)(f * (max-min));
}
bool bRand() {
    return (double)rand() / RAND_MAX < 0.5;
}

// Equals
bool test_case_equals(__float128 l, __float128 r, bool result) {
    bool a;
    a = l.Equals(r) == result;
    a = __float128::Equals(l, r) == result;
    return a;
}
void test_equals() {
    for(int i = 0; i < 10000; i++) {
        bool e = bRand();
        __float128 l = fRand(0.0, 1e15);
        __float128 r = l - !e * fRand(-1e30, 1e30);

        if (!test_case_equals(l, r, e)) {
            double left = 0.0, right = 0.0;
            __float128::ToDouble(l, left);
            __float128::ToDouble(r, right);
            std::cout << "fail   left:" << left << ", right:" << right << ", eql:" << e << std::endl;
        }
    }
}

// Compare
bool test_case_cmp(__float128 l, __float128 r, int result) {
    bool a;
    a = l.Cmp(r) == result;
    a = __float128::Cmp(l, r) == result;
    return a;
}
void test_cmp() {
    for(int i = 0; i < 10000; i++) {
        int c = iRand(-1, 1);
        __float128 l = fRand(0.0, 1e15);
        __float128 r = l - c * fRand(0, 1e30);
        
        if (!test_case_cmp(l, r, c)) {
            double left = 0.0, right = 0.0;
            __float128::ToDouble(l, left);
            __float128::ToDouble(r, right);
            std::cout << "fail   left:" << left << ", right:" << right << ", cmp:" << c << std::endl;
        }
    }
}

int main(char** argv, int argc) {
    test_equals();
    test_cmp();
    return 0;
}
