import os
import shutil
import glob
import pandas as pd
import re
from typing import List, Tuple
from pathlib import Path

def safe_copy_dir(src_dir: str, dst_dir: str, overwrite=False):
    """
    安全拷贝整个目录
    :param src_dir: 源目录路径 
    :param dst_dir: 目标路径
    :param overwrite: 是否合并覆盖（Python 3.8+ required）
    """
    src = Path(src_dir)
    dst = Path(dst_dir)
    
    if not dst.exists():
        shutil.copytree(src, dst)
    else:
        print(f"⏩ Skipping existing directory: {dst}")


def copypict(pattern,outputfile,templatedir="part04Template.md"):
    current_dir = os.getcwd()
    target_dir = os.path.join(current_dir, templatedir,"_v_images")
    source_file = glob.glob(pattern)[0]
    if os.path.exists(source_file):
        shutil.copy(source_file, os.path.join(target_dir,outputfile))
        print("Files successfully moved to the _v_images directory")
    else:
        print(f"Error: Source file {source_file} does not exist")

def update_csq(template_path: str, 
              new_table_data: List[Tuple[str, int, float]], 
              output_path: str = None) -> None:
    """
    精准更新Markdown报告中的CSQ表格，保持其他内容不变
    
    Args:
        template_path: 模板文件路径
        new_table_data: 新表格数据，格式为 [(label, count, frequency), ...]
        output_path: 输出路径（默认覆盖原文件）
    
    Returns:
        None
    """
    templatedir = os.path.dirname(template_path)
    save_to = os.path.join(templatedir, output_path) if output_path else template_path
    
    # 读取模板内容
    with open(template_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 1. 生成新的表格Markdown
    new_table_md = "\n".join([
        "|  lable   | count | frequency(%) |",
        "| :------: | :---: | :----------: |",
        *[f"| {row[0]:<8} | {row[1]:>5} | {row[2]:>12.1f} |" for row in new_table_data],
        ""  # 确保末尾有空行
    ])
    new_table_md =new_table_md+"\n"

    # 2. 精准定位表格区域（改进版正则）
    table_pattern = re.compile(
        r'(^|\n)(\|.*?lable.*?\|\n.*?)'  # 匹配表头
        r'(\|.*?\|.*?\|.*?\|.*?\n)+'     # 匹配表格内容行
        r'(\s*\n)?',                     # 匹配表格后的空行
        flags=re.MULTILINE
    )
    
    if not table_pattern.search(content):
        raise ValueError("No CSQ table found in the document. Check the template format")
    
    # 3. 精准替换（保留前后的内容）
    new_content = table_pattern.sub(r'\1' + new_table_md, content)
    
    # 4. 写入文件
    with open(save_to, 'w', encoding='utf-8') as f:
        f.write(new_content)
    
    print(f"✔ Table updated: {save_to}")

def generate_report_tables(CSQ2spdf, output_file=None):
    if not output_file:
        print("output_file was not specified")
        return

    # 按CSQ降序排序，取前10个高CSQ得分基因
    top_high = CSQ2spdf.sort_values('CSQ', ascending=False).head(10)
    
    # 按CSQ升序排序，取前10个低CSQ得分基因
    top_low = CSQ2spdf.sort_values('CSQ', ascending=True).head(10)
    
    def df_to_markdown(df):
        """将DataFrame转换为Markdown格式表格"""
        # 生成对齐方式行
        alignments = []
        for _ in df.columns:
            # 对GeneA和GeneB使用居中对齐，其他列右对齐
            if _ in ['GeneA', 'GeneB']:
                alignments.append(':---------:')
            else:
                alignments.append(':--------:')
        
        # 生成表头行
        header = '| ' + ' | '.join(df.columns) + ' |'
        # 生成对齐行
        alignment = '|' + '|'.join(alignments) + '|'
        # 生成数据行
        data_rows = []
        for _, row in df.iterrows():
            data_rows.append('| ' + ' | '.join(row.astype(str)) + ' |')
        
        return [header, alignment] + data_rows
    
    # 转换为带换行符的单一字符串
    high_table_str = '\n'.join(df_to_markdown(top_high))
    low_table_str = '\n'.join(df_to_markdown(top_low))
    
    # ------------------------------
    # 核心修改区：读取文件 -> 正则替换表格 -> 覆盖写回
    # ------------------------------
    try:
        with open(output_file, 'r', encoding='utf-8') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"❌ File not found: {output_file}")
        return

    # 定义你的定位说明文字（确保和 Markdown 模板里的文字完全一致）
    high_header = "列取CSQ得分最高的10个同源基因对情况（即同源基因具有相似的剪切模式）："
    low_header = "列取CSQ得分最低的10个同源基因对情况（即同源基因具有较大差异的剪切模式）："

    # 编译正则表达式
    # 逻辑：匹配 "你的定位文字 + 换行" (\g<1>)，然后匹配 "紧跟其后的所有以 | 开头的表格行"
    # re.MULTILINE 确保 ^ 能匹配每一行的开头
    pattern_high = re.compile(r'(' + re.escape(high_header) + r'\s*\n+)(?:^[ \t]*\|.*$\n?)*', re.MULTILINE)
    pattern_low = re.compile(r'(' + re.escape(low_header) + r'\s*\n+)(?:^[ \t]*\|.*$\n?)*', re.MULTILINE)

    # 进行替换：保留原标题（\g<1>），然后拼接你新生成的高分/低分表格，并加上换行
    content = pattern_high.sub(r'\g<1>' + high_table_str + '\n\n', content)
    content = pattern_low.sub(r'\g<1>' + low_table_str + '\n\n', content)

    # 重新将更新后的完整内容覆盖写回文件
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(content)

    print(f"✔ Table updated: {output_file}")

def count_lines(filepath):
    """快速计数文件行数（不加载内容）"""
    with open(filepath, 'r', encoding='utf-8') as f:
        return sum(1 for _ in f)

def replaceKeywordNumber(num,keyward, template_path):
    # 1. 读取原文件全部内容
    with open(template_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # 2. 编译正则表达式
    pattern = re.compile(rf'(\*\s*{keyward}[：:]\s*)\d+')

    if pattern.search(content):
        # 3. 替换数字：保留匹配到的前缀 \g<1>，把后面的 \d+ 替换为传入的 num
        new_content = pattern.sub(rf'\g<1>{num}', content)
        
        # 4. 覆盖写回原文件
        with open(template_path, 'w', encoding='utf-8') as f:
            f.write(new_content)
        
        print(f"✔ Successfully updated {keyward} count to: {num}")
    else:
        print(f"⚠️ Warning: In {template_path}, no text matching '* {keyward}:[number]' was found; no changes made.")


def replace_markdown_table(table_data, keyword, template_path):
    import re
    
    #print(f"✔ 正在处理关键字: {keyword}")
    #print(f'替换内容为:{table_data}')
    # 正则表达式匹配表格中的具体行（确保只匹配指定行）
    table_pattern = re.compile(
        rf'^(\|\s*{re.escape(keyword)}\s*\|)\s*(.+?)\s*(\|.*)$',  # 捕获目标行
        flags=re.MULTILINE
    )


    # 读取模板 Markdown 内容
    with open(template_path, "r", encoding="utf-8") as f:
        content = f.read()

    # 测试打印文件内容
    #print("✔ 模板文件内容:")
    #print(content)

    # 检查是否匹配到指定的表格内容
    match = table_pattern.search(content)
    if not match:
        raise ValueError(f"❌ No table section found for keyword '{keyword}'. Check the template file contents")

    #print(f"✔ 找到匹配的行: {match.group()}")  # 打印找到的部分

    # 替换匹配的表格行中的数字和百分比
    updated_content = table_pattern.sub(
        rf'\1 {table_data} \3',  # 保留组 1 和组 3，替换组 2 的内容
        content
    )


    # 输出文件路径（这里覆盖原模板文件）
    save_path = template_path
    with open(save_path, "w", encoding="utf-8") as f:
        f.write(updated_content)

    #print(f"✔ 表格内容已替换并保存到: {save_path}")

def merge_report():
    # 获取当前工作目录
    base_path = Path.cwd()
    
    # 定义输出目录和文件
    out_dir = base_path / "merge_report"
    out_img = out_dir / "_v_images"
    out_md = out_dir / "merged_report.md"
    
    # 创建输出目录和图片子目录（parents=True会自动创建多级目录，exist_ok=True表示存在不报错）
    out_img.mkdir(parents=True, exist_ok=True)
    
    # 规定合并的目录顺序
    sub_dirs = ["CSQ", "CSQ2sp", "RNAseq"]
    
    # 使用写模式 "w" 打开目标 md 文件（这会清空之前可能存在的旧文件）
    with open(out_md, "w", encoding="utf-8") as f_out:
        for folder in sub_dirs:
            # 【修改点 1】: target_dir 应该是目录，而不是直接拼上 "Template.md"
            target_dir = base_path / folder  / "Template.md"
            
            # 判断目标目录是否存在
            if target_dir.is_dir():
                print(f"Processing directory: {target_dir}")
                
                # 1. 查找并合并所有的 .md 文件 (这里会自动找到 Template.md 等文件)
                for md_file in target_dir.glob("*.md"):
                    if md_file.is_file():
                        print(f"  -> Merging file: {md_file.name}")
                        with open(md_file, "r", encoding="utf-8") as f_in:
                            content = f_in.read()
                            f_out.write(content)
                            # 在每个文件末尾追加两个换行符，防止内容紧贴着粘连
                            f_out.write("\n\n") 
                            
                # 2. 拷贝 _v_images 目录内的所有图片
                img_dir = target_dir / "_v_images"
                if img_dir.is_dir():
                    print(f"  -> Copying images: from {img_dir} to {out_img}")
                    for img_file in img_dir.iterdir():
                        if img_file.is_file(): 
                            shutil.copy2(img_file, out_img)
                else:
                    print(f"  -> Image directory not found: {img_dir}") # 增加一条日志，方便调试
            else:
                print(f"Skipping: {target_dir} does not exist")
    print("-" * 40)
    print("Merge complete!")
    print(f"Final report directory: {out_dir}")
    print(f"Merged Markdown: {out_md}")

# 测试运行
def main(args_list):
    merge_report()
    

if __name__ == '__main__':
    main()
    
