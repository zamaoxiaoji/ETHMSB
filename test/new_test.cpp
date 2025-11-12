#include <iostream>
#include <cmath>
#include <cstdint>
#include <bitset> 

int main() {
    const uint16_t pow5_max = (1 << 5);   
    const uint16_t pow10_max = (1 << 10); 
    const uint16_t pow11_max = (1 << 11); 
    const uint32_t div_factor = (1 << 16); 

    for (uint16_t high5 = 0; high5 < pow5_max; ++high5) {
        uint16_t m = (high5 << 11); // 高 8 位，低 8 位为 0
        uint16_t ectract = (m & (1 << 12)); 

        for (int16_t e = -(pow10_max); e < pow10_max; ++e) {
            uint16_t m1 = m  + e + (1<<10);

            // m1c / 2^16 并四舍五入
            // uint16_t m2 = static_cast<uint16_t>(std::round(static_cast<double>(m1) / div_factor));
            // uint16_t m3 = static_cast<uint16_t>(std::round(static_cast<double>(m) / div_factor));

            uint16_t m2 = (m1 & (1 << 15))/(1 << 15); 
            uint16_t m3 = (m & (1 << 15))/(1 << 15);
            
            auto print16 = [](uint32_t v) {           // 捕获值，统一转成 16 位
                std::cout << std::bitset<16>(v) << '\n';
            };

            
            if (m2 != m3) {
                std::cout << "Mismatch found:\n";
                std::cout << "m = " << m << ", e = " << e << ", m1 = " << m1 << "\n";
                std::cout << "m : ";
                print16(m);       // 打印 m 的 16 位二进制
                std::cout << "m-: ";
                print16(m - ectract);
                std::cout << "e : ";
                print16(e);       // 打印 e 的 16 位二进制（高位用 0 补齐）
                print16(m1);     // 打印 m1c 的低 16 位（二进制序列）
                
                std::cout << "m2 = " << m2 << ", m3 = " << m3 << "\n\n";
            }
        }
    }

    std::cout << "ok\n";
    return 0;
}
