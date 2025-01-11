from concurrent.futures import ProcessPoolExecutor
import os
import subprocess
import re
from tqdm import tqdm
import sys
from multiprocessing import Manager

output_folder = "txt_folder"
os.makedirs(output_folder, exist_ok=True)

ell = int(sys.argv[1])
success = 10
problem_id = 6

def run_command(args):
    ell, index, results = args  # 解包參數
    command = f"./sweep {ell} {success} {problem_id} {index}"
    
    try:
        # 執行指令並捕獲輸出
        output = subprocess.check_output(command, shell=True, text=True)

        # 使用正則表達式提取 NFE
        match = re.search(r"NFE:\s*([\d.]+)", output)
        if match:
            nfe = float(match.group(1))
            results[index] = nfe  # 更新共享字典

            # 儲存結果到對應的檔案
            output_file = os.path.join(output_folder, f"{ell}_{index}.txt")
            with open(output_file, "w") as f:
                f.write(f"ell: {ell}, NFE: {nfe}\n")
        else:
            print(f"未找到 NFE 值: {output}")
    except subprocess.CalledProcessError as e:
        print(f"執行失敗: {e}")

# 確保 ell 和其他參數事先設置好
count = 1000

# 使用 Manager 共享結果字典
with Manager() as manager:
    results = manager.dict()  # 使用共享字典
    tasks = [(ell, index, results) for index in range(1, count + 1)]
    
    # 使用 ProcessPoolExecutor 進行平行處理
    with ProcessPoolExecutor() as executor:
        list(tqdm(executor.map(run_command, tasks), total=len(tasks)))

    # 計算平均 NFE（在所有任務完成後）
    if len(results) == count:  # 確保所有任務都完成
        average_nfe = sum(results.values()) / len(results) if results else 0
        print(f"ell: {ell}, 平均NFE: {average_nfe}")
    else:
        print(f"部分任務未完成，完成數量: {len(results)}/{count}")
