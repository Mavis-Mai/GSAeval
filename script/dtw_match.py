import numpy as np
import heapq
from collections import defaultdict

def phase_matrix(x, y,window=5):
    n, m = len(x), len(y)
    if n == 0 or m == 0:
        return []
    
    # 1. 计算成本矩阵（处理None值）
    cost_matrix = np.zeros((n, m))
    for i in range(n):
        for j in range(max(0, i-window), min(m, i+window+1)):
            cost_matrix[i,j] = np.inf if (x[i] is None or y[j] is None) \
            else abs(x[i]-y[j])
    return cost_matrix

def len_matrix(x, y, window=5, max_cost=1.5, max_skip=1):
    """
    Args:
        x, y: 输入序列
        window: 匹配索引的最大偏移
        max_cost: 允许的单点最大成本（超过则设为inf）
        max_skip: 最大连续跳过次数（避免跳过过多点）
    """
    n, m = len(x), len(y)
    if n == 0 or m == 0:
        return np.zeros((n, m))
    
    cost_matrix = np.full((n, m), np.inf)
    for i in range(n):
        # 动态调整窗口大小（如果附近有异常值）
        current_window = window
        if i > 0 and cost_matrix[i-1].min() == np.inf:
            current_window = min(window + 2, m)  # 扩大窗口
        
        for j in range(max(0, i-current_window), min(m, i+current_window+1)):
            if x[i] is None or y[j] is None:
                continue
                
            diff = abs(x[i] - y[j])
            denominator = (x[i] + y[j]) / 2
            
            if denominator == 0:
                cost = diff
            else:
                cost = diff / denominator
            
            # 如果成本超过阈值，则跳过（设为inf）
            if cost > max_cost:
                continue
            cost_matrix[i, j] = cost
    
    return cost_matrix

def find_unique_top_matches(cost_matrix,window=5, top_k=3):
    """优化版本：保证唯一匹配路径且允许缺失匹配"""
    n,m=cost_matrix.shape
    # 2. 扩展的动态规划，包含缺失匹配选项
    dp = np.full((n+1, m+1), np.inf)
    dp[0,0] = 0
    
    # 记录路径来源，用于回溯
    trace = defaultdict(list)
    trace[(0,0)] = []
    
    for i in range(n+1):
        for j in range(m+1):
            if i == 0 and j == 0:
                continue
                
            options = []
            # 匹配选项（斜对角）
            if i > 0 and j > 0 and abs(i-j) <= window:
                match_cost = dp[i-1,j-1] + cost_matrix[i-1,j-1]
                options.append((match_cost, (i-1,j-1), True))  # 最后布尔值表示是否匹配
            
            # x序列缺失选项（向上）
            if i > 0:
                skip_cost = dp[i-1,j] + 1  # 固定缺失成本
                options.append((skip_cost, (i-1,j), False))
                
            # y序列缺失选项（向左）
            if j > 0:
                skip_cost = dp[i,j-1] + 1  # 固定缺失成本
                options.append((skip_cost, (i,j-1), False))
            
            if options:
                min_cost, prev, is_match = min(options, key=lambda x: x[0])
                dp[i,j] = min_cost
                trace[(i,j)] = trace[prev].copy()
                if is_match:
                    trace[(i,j)].append((i-1, j-1))  # 只记录实际匹配点
    
    # 3. 从终点回溯所有可能路径
    unique_paths = set()
    
    def get_unique_path_key(path):
        """生成路径的唯一标识，忽略顺序"""
        return frozenset(path)  # 使用集合来避免顺序影响
    
    if (n,m) in trace:
        path = trace[(n,m)]
        unique_paths.add(get_unique_path_key(path))
    
    # 4. 收集top_k条最佳唯一路径
    result = []
    for path_key in sorted(unique_paths, key=lambda k: sum(cost_matrix[i,j] for i,j in k)):
        if len(result) >= top_k:
            break
        # 恢复原始顺序（按x索引排序）
        sorted_path = sorted(path_key, key=lambda x: x[0])
        result.append(sorted_path)
    
    return result


