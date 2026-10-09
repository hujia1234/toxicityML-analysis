import os
import sys
import glob

def find_and_rename_files():
    print("🔍 正在当前文件夹中寻找数据文件...")

    gse_files = glob.glob("GSE*.txt")
    if len(gse_files) == 0:
        print("❌ 错误：未找到以 GSE 开头的 .txt 文件！")
        return False
    elif len(gse_files) > 1:
        print(f"⚠️  警告：找到多个 GSE 开头的文件：{gse_files}")
        print(f"    将使用第一个：{gse_files[0]}")

    gse_file = gse_files[0]

    gpl_files = glob.glob("GPL*.txt")
    if len(gpl_files) == 0:
        print("❌ 错误：未找到以 GPL 开头的 .txt 文件！")
        return False
    elif len(gpl_files) > 1:
        print(f"⚠️  警告：找到多个 GPL 开头的文件：{gpl_files}")
        print(f"    将使用第一个：{gpl_files[0]}")

    gpl_file = gpl_files[0]

    try:
        if os.path.exists("probe.txt"):
            os.remove("probe.txt")
        if os.path.exists("ann.txt"):
            os.remove("ann.txt")

        os.rename(gse_file, "probe.txt")
        os.rename(gpl_file, "ann.txt")
        print(f"✅ 文件重命名成功：")
        print(f"   {gse_file} -> probe.txt")
        print(f"   {gpl_file} -> ann.txt")
        return True
    except Exception as e:
        print(f"❌ 重命名失败：{str(e)}")
        return False

def find_gene_symbol_column(header_line):
    headers = header_line.split('\t')
    clean_headers = [h.replace('"', '').replace("'", '').strip() for h in headers]

    print(f"\n✅ 识别到的表头列名：{clean_headers}")

    match_rules = [
        lambda h: 'gene' in h.lower() and 'sym' in h.lower(),
        lambda h: h.lower() in ['symbol', 'sym', 'gene_symbol', 'genesymbol', 'gene_sym'],
        lambda h: 'gene' in h.lower() and 'name' in h.lower()
    ]

    for rule_idx, rule in enumerate(match_rules):
        for col_idx, header in enumerate(clean_headers):
            if rule(header):
                print(f"\n✅ 自动识别到 Gene Symbol 列：第 {col_idx+1} 列，列名：{clean_headers[col_idx]}")
                return col_idx
    return -1

def main():
    if not find_and_rename_files():
        return

    exp_file = "probe.txt"
    gpl_file = "ann.txt"
    output_file = "geneMatrix.txt"

    gene_symbol_col = -1
    probe_gene_hash = {}
    header_found = False

    print(f"\n正在读取 {gpl_file}，跳过注释行，寻找表头...")

    with open(gpl_file, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('!'):
                continue

            if not header_found:
                print(f"\n✅ 在第 {line_num} 行找到表头")
                gene_symbol_col = find_gene_symbol_column(line)

                if gene_symbol_col == -1:
                    print("\n⚠️  自动识别失败，请手动输入Gene Symbol所在的列号（从1开始数）：")
                    try:
                        gene_symbol_col = int(input().strip()) - 1
                    except ValueError:
                        print("❌ 输入错误，请输入有效的数字！")
                        return

                header_found = True
                continue

            gpl_arr = line.split('\t')
            if (gene_symbol_col < len(gpl_arr) and
                gpl_arr[gene_symbol_col] != '' and
                ' ' not in gpl_arr[gene_symbol_col] and
                not gpl_arr[gene_symbol_col].startswith('---') and
                not gpl_arr[gene_symbol_col][0].isdigit()):

                gene_name = gpl_arr[gene_symbol_col]
                if '///' in gene_name:
                    gene_name = gene_name.split('///')[0].strip()
                gene_name = gene_name.replace('"', '').replace("'", '')

                probe_id = gpl_arr[0].strip()
                probe_gene_hash[probe_id] = gene_name

    if len(probe_gene_hash) == 0:
        print("❌ 错误：没有读取到有效的探针-基因对应关系")
        return
    print(f"\n✅ 成功读取到 {len(probe_gene_hash)} 个探针的基因注释")

    print(f"\n正在读取 {exp_file} 表达矩阵...")
    sample_names = []
    data_hash = {}

    with open(exp_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('!'):
                continue

            line = line.replace('"', '')
            exp_arr = line.split('\t')

            if line.startswith('ID_REF'):
                sample_names.append("ID_REF")
                for i in range(1, len(exp_arr)):
                    raw_sample_name = exp_arr[i]
                    if '_' in raw_sample_name:
                        clean_name = raw_sample_name.split('_')[0]
                    elif '.' in raw_sample_name:
                        clean_name = raw_sample_name.split('.')[0]
                    else:
                        clean_name = raw_sample_name
                    sample_names.append(clean_name)
                print(f"✅ 识别到 {len(sample_names)-1} 个样本")
            else:
                for i in range(len(sample_names)):
                    val = exp_arr[i].replace('"', '').strip() if i < len(exp_arr) else ""
                    if sample_names[i] not in data_hash:
                        data_hash[sample_names[i]] = []
                    data_hash[sample_names[i]].append(val)

    if "ID_REF" not in data_hash:
        print("❌ 错误：probe.txt格式不正确")
        return
    probe_names = data_hash["ID_REF"]
    del data_hash["ID_REF"]

    print("\n正在计算基因平均表达量...")
    gene_list = set()
    sample_gene_exp = {}

    for sample_name in data_hash:
        exp_values = data_hash[sample_name]
        gene_count = {}
        gene_sum = {}

        for i in range(len(probe_names)):
            probe_id = probe_names[i].strip()
            if probe_id in probe_gene_hash:
                gene_name = probe_gene_hash[probe_id]
                gene_list.add(gene_name)

                if gene_name not in gene_count:
                    gene_count[gene_name] = 0
                    gene_sum[gene_name] = 0.0

                try:
                    val = float(exp_values[i])
                except (ValueError, IndexError):
                    val = 0.0

                gene_count[gene_name] += 1
                gene_sum[gene_name] += val

        gene_average = {}
        for gene in gene_count:
            gene_average[gene] = gene_sum[gene] / gene_count[gene]
        sample_gene_exp[sample_name] = gene_average

    sample_names[0] = "geneNames"
    sorted_genes = sorted(gene_list)

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write('\t'.join(sample_names) + '\n')
        for gene in sorted_genes:
            line_content = [gene]
            for i in range(1, len(sample_names)):
                sample = sample_names[i]
                exp_val = sample_gene_exp[sample].get(gene, "")
                line_content.append(str(exp_val))
            f.write('\t'.join(line_content) + '\n')

    print(f"\n🎉 全部处理完成！")
    print(f"✅ 最终生成基因数：{len(sorted_genes)}")
    print(f"✅ 结果文件已保存为：{output_file}")

if __name__ == "__main__":
    main()
