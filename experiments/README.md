# experiments 参数说明

本目录下的实验程序都属于“只新增文件、复制函数到新命名空间插桩”的方式实现，不修改 HE3DB/TFHEpp 既有源码。

## 构建方式（CMake）

在项目根目录执行：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

生成的可执行文件默认在 `build/bin/`。

注意：实验程序为了“失败样本同 seed 复跑”会调用 `TFHEpp::generator.seed(seed)`，因此需要 `USE_RANDEN=ON`（本仓库 TFHEpp 默认开启）。

---

## 1) `msb5_fail_trace`

用途：验证对 `plain_bits > 5` 的多值整数密文，直接套用 PBS-based 的 `ExtractMSB5` 会在 worst-case（接近回绕边界、长串 1）输入上出现错误；并在失败样本复跑时打印三段格式 phase 对比：

```
<m区> | <Δ间隔区> | <误差区>
```

可执行文件：`build/bin/msb5_fail_trace`

### 参数列表

- `--trials N`：每个 `plain_bits` 的试验次数（默认 `500`）
- `--min bits`：扫描的 `plain_bits` 起始值（默认 `6`）
- `--max bits`：扫描的 `plain_bits` 结束值（默认 `10`）
- `--err-bits e`：三段格式里“误差区”的低位 bit 数（默认 `8`；可设 `3` 复现你给的示例风格）
- `--case msb0|msb1`：worst-case 明文选择（默认 `msb1`）
  - `msb0`：`m = 2^(t-1) - 1`（形如 `0111..11`）
  - `msb1`：`m = 2^t - 1`（形如 `1111..11`，最接近回绕边界）
- `--K K`：每个 `plain_bits` 最多打印 K 个失败样本 trace（默认 `5`）
- `--base-seed S`：随机种子基数（支持十进制或 `0x...`，默认 `0xC0FFEE`）
- `--no-control`：不跑对照输入 `m3 = 3*2^(t-2)`（默认会跑）
- `--stop`：某个 `plain_bits` 出现失败后停止继续扫描
- `--help` / `-h`：打印帮助

### 运行示例

```bash
./build/bin/msb5_fail_trace --trials 2000 --min 10 --max 10 --K 1 --stop --err-bits 8 --case msb1
```

---

## 2) `msb_wraparound_experiment`

用途：TFHEpp 层面的 PBS/MSB wrap-around 实验（使用 `experiments/tfhepp_instrumented_pbs.hpp` 的复制版接口）。

可执行文件：`build/bin/msb_wraparound_experiment`

### 参数列表

- `--trials N`：每个 t 的试验次数（默认 `2000`）
- `--t-min N`：扫描的 t 起始值（默认 `8`）
- `--t-max N`：扫描的 t 结束值（默认 `32`）
- `--alpha-in X`：输入加密噪声（默认 `0.0`）
- `--fail-print-limit K`：每个 t 最多打印 K 个失败样本（默认 `5`）
- `--no-stop`：不要在首次出现失败的 t 停止
- `--no-control`：跳过对照输入（默认会跑）
- `--help` / `-h`：打印帮助

---

## 3) `msb6_test`

用途：对比 **原始** `ExtractMSB5` 与 **实验版** `MSB6(offset=Δ/2)` 在 **6-bit 输入**（默认 `m=63(111111)`）上的失败率。

可执行文件：`build/bin/msb6_test`

### 参数列表

- `--trials N`：试验次数（默认 `500`）
- `--base-seed S`：随机种子基值（默认 `0xC0FFEE`）
- `--case msb0|msb1`：默认 `msb1`
  - `msb1`：`m=63(111111)`（最坏情况，最接近回绕边界）
  - `msb0`：`m=31(011111)`
- `--m M`：指定明文 `m`（覆盖 `--case`）

### 运行示例

```bash
./build/bin/msb6_test --trials 500 --case msb1
```

---

## 4) `msb_timing_bench`

用途：对不同 **合法** `plain_bits` 的 MSB 提取路径（`HomMSB` 内部会分派到 `ExtractMSB5/9/10` 或 `ImExtractMSB*`）做运行时间对比，输出每次调用的耗时统计（平均/中位数/min/max）。

可执行文件：`build/bin/msb_timing_bench`

### 参数列表

- `--level 1|2`：输入密文等级
  - `1`：输入为 `TLWELvl1`（合法 `plain_bits<=10`）
  - `2`：输入为 `TLWELvl2`（合法 `plain_bits<=33`）
- `--min bits` / `--max bits`：扫描的 `plain_bits` 范围
- `--iters N`：每个 `plain_bits` 的计时迭代次数（默认 `20`）
- `--warmup N`：预热次数（默认 `2`，不计时）
- `--base-seed S`：用于生成输入密文的随机种子基数（默认 `0xC0FFEE`）
- `--also-msb6`：在 `level=1 & plain_bits=6` 时，额外测一次实验版 `MSB6(offset=Δ/2)`（便于和 `HomMSB(ExtractMSB9)` 对比）
- `--csv`：输出 CSV（便于画图）

### 运行示例

```bash
# lvl1: plain_bits 5..10
./build/bin/msb_timing_bench --level 1 --min 5 --max 10 --iters 20 --warmup 2

# lvl1: 额外比较 MSB6
./build/bin/msb_timing_bench --level 1 --min 6 --max 6 --iters 50 --also-msb6

# lvl2: plain_bits 5..33（注意更慢）
./build/bin/msb_timing_bench --level 2 --min 5 --max 33 --iters 10 --warmup 1

# 输出 CSV
./build/bin/msb_timing_bench --level 1 --min 5 --max 10 --iters 20 --csv
```
