import subprocess
import re
import os
from concurrent.futures import ProcessPoolExecutor
import sys

continuousSuccess = sys.argv[1] # 連續幾次成功
problem = sys.argv[2]
ell = [36, 72, 108, 144]

print(f"problem id = {problem}")
# 0 onemax
# 1 mk
# 2 FTRAP
# 3 cyc
# 4 nk
# 5 SPIN
# 6 SAT
# 7 mixTrap
# 8 3-opt

# 結果存放資料夾
output_folder = "txt_folder"
os.makedirs(output_folder, exist_ok=True)

# 存放結果的字典
results = {}
def run_command(ell):
    """執行指令，提取 NFE，並將結果存檔"""
    if problem == 6:
        command = f"./sweep {ell} {continuousSuccess} {problem} 0"
    else:
        command = f"./sweep {ell} {continuousSuccess} {problem}"

    try:
        # 執行指令並捕獲輸出
        output = subprocess.check_output(command, shell=True, text=True)

        # 使用正則表達式提取 NFE
        match = re.search(r"NFE:\s*([\d.]+)", output)
        if match:
            nfe = float(match.group(1))
            results[ell] = nfe

            # 儲存結果到對應的檔案
            output_file = os.path.join(output_folder, f"{ell}.txt")
            with open(output_file, "w") as f:
                f.write(f"ell: {ell}, NFE: {nfe}\n")

            print(f"ell: {ell}, NFE: {nfe}")
        else:
            print(f"未找到 NFE 值: {output}")
    except subprocess.CalledProcessError as e:
        print(f"執行失敗: {e}")

# 使用 ProcessPoolExecutor 進行平行執行
with ProcessPoolExecutor() as executor:
    executor.map(run_command, ell)
