#include <iostream>
#include <cmath>
#include <cstdint>
#include <bitset> 
#include <iomanip>

int main() {
    const uint16_t pow9_max = (1 << 9);   
    const uint16_t pow10_max = (1 << 10); 
    const uint16_t pow11_max = (1 << 11); 
    const uint32_t div_factor = (1 << 16); 

    // 一共16bit，高9比特是消息，考虑原始msb能够提取5bit的。
    for (uint16_t high9 = 0; high9 < pow9_max; ++high9) {
        uint16_t m = (high9 << 7);

        uint16_t ectract = (m & (1 << 11)); 
        // 
        for (int16_t e = -(pow10_max); e < pow10_max; ++e) {
            uint16_t m1 = m - ectract + e + (1<<10);

            // 除 2^16 并四舍五入
            // uint16_t m2 = static_cast<uint16_t>(std::round(static_cast<double>(m1) / div_factor));
            // uint16_t m3 = static_cast<uint16_t>(std::round(static_cast<double>(m) / div_factor));
            
            uint16_t m2 = (m1 & (1 << 15))/(1 << 15); 
            uint16_t m3 = (m & (1 << 15))/(1 << 15);

            auto print16 = [](uint32_t v) {    // 转16位二进制输出
                std::cout << std::bitset<16>(v) << '\n';
            };

            
            if (m2 != m3) {
                std::cout << "Mismatch found:\n";
                std::cout << "m = " << m << ", e = " << e << ", m1 = " << m1 << "\n";
                std::cout << "m : ";
                print16(m);       
                std::cout << "m-: ";
                print16(m - ectract);
                std::cout << "e : ";
                print16(e);  
                std::cout << "m1 : ";     
                print16(m1);     
                
                std::cout << "m2 = " << m2 << ", m3 = " << m3 << "\n\n";
            }
        }
    }
    std::cout << "ok \n" ;

    return 0;

}
