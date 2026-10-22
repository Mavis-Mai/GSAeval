import numpy as np

class jaccard_tools:
    def get_kmer_set(self, sequence, k=3):
        """将序列切割为k-mer集合"""
        return {sequence[i:i+k] for i in range(len(sequence)-k+1)}

    def jaccard_similarity(self, seq1, seq2, k=3):
        """计算Jaccard相似度"""
        set1 = self.get_kmer_set(seq1, k)
        set2 = self.get_kmer_set(seq2, k)
        intersection = set1 & set2
        union = set1 | set2
        return len(intersection) / len(union) if union else 0

    def build_jaccard_matrix(self, list_A, list_B, k=3):
        """构建全连接的Jaccard相似度矩阵"""
        m, n = len(list_A), len(list_B)
        matrix = np.zeros((m, n))
        for i in range(m):
            for j in range(n):
                matrix[i][j] = self.jaccard_similarity(list_A[i], list_B[j], k)
        return matrix

    def find_constrained_path_fixed(self, jaccard_matrix):
        m, n = jaccard_matrix.shape
        dp = -np.inf * np.ones((m, n))
        path_pre = -1 * np.ones((m, n), dtype=int)
        valid_rows = []  # 记录有效行索引（非全0行）

        # 初始化：标记有效行
        for i in range(m):
            if np.any(jaccard_matrix[i] > 0):
                valid_rows.append(i)

        if not valid_rows:
            return [], 0.0  # 全矩阵为0

        # 初始化第一有效行
        first_row = valid_rows[0]
        for j in range(n):
            if jaccard_matrix[first_row, j] > 0:
                dp[first_row, j] = jaccard_matrix[first_row, j]

        # 动态规划：仅处理有效行
        for k in range(1, len(valid_rows)):
            i = valid_rows[k]      # 当前有效行
            prev_i = valid_rows[k-1]  # 前一有效行
            for j in range(n):
                if jaccard_matrix[i, j] == 0:
                    continue
                for prev_j in range(j):
                    if dp[prev_i, prev_j] + jaccard_matrix[i, j] > dp[i, j]:
                        dp[i, j] = dp[prev_i, prev_j] + jaccard_matrix[i, j]
                        path_pre[i, j] = prev_j

        # 回溯路径
        max_score = -np.inf
        last_row, last_col = -1, -1
        for i in valid_rows:
            for j in range(n):
                if dp[i, j] > max_score:
                    max_score = dp[i, j]
                    last_row, last_col = i, j

        if max_score == -np.inf:
            return [], 0.0

        # 反向追踪路径
        path = []
        current_row, current_col = last_row, last_col
        while current_row >= 0 and current_col >= 0:
            path.append((current_row, current_col))
            current_col = path_pre[current_row, current_col]
            current_row = valid_rows[valid_rows.index(current_row)-1] if current_row != valid_rows[0] else -1
        path.reverse()

        return path, max_score