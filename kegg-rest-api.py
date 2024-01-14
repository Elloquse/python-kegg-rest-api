from matplotlib import pyplot as plt
from scipy.stats.stats import pearsonr
from scipy.stats import shapiro

class Gene():
    
    def __init__(self, entry, organism, motif, aaseq, ntseq):
        self.entry = entry
        self.organism = organism
        self.motif = motif
        self.aaseq = aaseq
        self.ntseq = ntseq
        
        
    def kegg_gene(gene_id):
        import requests
        import re
        
        adress = f'http://rest.kegg.jp/get/{gene_id}'
        resp = requests.get(url=adress)
        all_lines = resp.text.splitlines()
         
        gene = [[],[],[],[],[]]
        for i, j in enumerate(all_lines):
            if 'ENTRY' in j:
                gene[0].append(re.split(r'\W+', j)[1])
            if 'ORGANISM' in j:
                gene[1].append(re.split(r'\W+', j)[1])
            if 'MOTIF' in j:
                gene[2].append(re.split(r': ', j)[1])
            if 'AASEQ' in j:
                aaseq_ind = i
            if 'NTSEQ' in j:
                ntseq_ind = i
            if '///' in j:
                end_ind = i
        aaseq_before = (''.join([seq for seq in all_lines[aaseq_ind+1:ntseq_ind]])).split()
        aaseq = ''
        for i in aaseq_before:
            aaseq += i
        gene[3].append(aaseq)
        
        ntseq_before = (''.join([seq for seq in all_lines[ntseq_ind+1:end_ind]])).split()
        ntseq = ''
        for i in ntseq_before:
            ntseq += i
        gene[4].append(ntseq)
        
        
        # Возврат пустых значений
        for i in gene:
            if len(i) == 0:
                i.append(None)
                
        return Gene(gene[0], gene[1], gene[2], gene[3], gene[4])
        

    
    def gene_ort(gene_id, gene_number):
        import requests
        import re
        
        adress = f'http://www.kegg.jp/ssdb-bin/ssdb_best?org_gene={gene_id}'
        resp = requests.get(url=adress)
        all_lines = resp.text.splitlines()
        
        for i, j in enumerate(all_lines):
            if '-----' in j:
                start_ind = i 
                break
        
        data_list = []
        for i in all_lines[start_ind+1:start_ind+gene_number+1]:
            if 'INPUT TYPE' in i:
                data_list.append(i)
            else:
                break
        

        ort = []
        for i in data_list:
            ort_gene = []
            ort_gene.append(re.search('(?<=TARGET=\"_blank\">)([^<\/>]+)', i).group())
            ort_gene.append((re.search('(?<=[0-9\s+])\d+(?=\s+\()', i).group()).split()[0])
            ort_gene.append((re.search('(?<=[0-9]\))([^-&]+)', i).group()).split()[1])
            ort.append(ort_gene)
        return ort
    
    
    def ort_motif(gene_id, gene_number):
        import requests
        import re
        from matplotlib import pyplot as plt
         
        url = Gene.gene_ort(gene_id, gene_number)
        
        motifs = []
        for i in url:
            adress = f'https://www.kegg.jp/ssdb-bin/ssdb_motif?kid={i[0]}'
            resp = requests.get(url=adress)
            all_lines = resp.text.splitlines()
            
            table_range = [i for i, j in enumerate(all_lines) if '</table>' in j]
            mots = []
            for j in all_lines[table_range[0]:table_range[1]]:
                if 'https://www.genome' in j:
                    mots.append(re.search('(?<=pf:)([^\"]+)', j).group())
            motifs.append(mots)

            
        adress = f'https://www.kegg.jp/ssdb-bin/ssdb_motif?kid={gene_id}'
        resp = requests.get(url=adress)
        all_lines = resp.text.splitlines()
        table_range = [i for i, j in enumerate(all_lines) if '</table>' in j]
        gene_motif = []
        for j in all_lines[table_range[0]:table_range[1]]:
            if 'https://www.genome' in j:
                gene_motif.append(re.search('(?<=pf:)([^\"]+)', j).group())
        motifs.append(gene_motif)

        
        # выявление всех существующих мотивов
        all_motifs = []
        for i in motifs: 
            for j in i:
                    all_motifs.append(j)
        all_motifs = list(set(all_motifs))

        # составление словаря
        dic_motif = {}
                       
        count = 0
        for i in all_motifs:
            for j in motifs:
                if i in j:
                    count += 1
            dic_motif[i] = count
            count = 0
        
        min_lim = min(dic_motif.values()) - (sum(dic_motif.values()) * 0.01)
        max_lim = max(dic_motif.values()) + (sum(dic_motif.values()) * 0.01)

        fig = plt.figure()
        ax = fig.add_axes([0,0,1,1])
        ax.bar(dic_motif.keys(), dic_motif.values(), color='#DF6B59')
        ax.set_ylim(min_lim, max_lim)
        ax.set_ylabel('Встречаемость мотивов', fontsize=14)
        plt.xticks(rotation=340, fontsize=10)

        return dic_motif
        
        

#your_gene = Gene.kegg_gene('hsa:5555')
#ort_100 = Gene.gene_ort('hsa:5555', 100)
#ort_10 = Gene.gene_ort('hsa:5555', 10)
#mot_100 = Gene.ort_motif_100('hsa:7314')
#mot_10 = Gene.ort_motif('hsa:5555', 10)

#### Пирсон ####

# ---- для ортологов ----

# x_100 = [float(i[1]) for i in ort_100]
# y_100 = [float(i[2]) for i in ort_100]
# p_100 = pearsonr(x_100, y_100)
# print(f'Коэффициент корреляции Пирсона для 100 ортологов: {p_100[0]}')
# print(f'Уровень значимости: {p_100[1]}')


# x_10 = [float(i[1]) for i in ort_10]
# y_10 = [float(i[1]) for i in ort_10]
# p_10 = pearsonr(x_10, y_10)
# print(f'Коэффициент корреляции Пирсона для 10 ортологов: {p_10[0]}')
# print(f'Уровень значимости: {p_10[1]}')


#### Шапиро-Уилк #####

# ---- для мотивов ----

# x_m_100 = [float(i) for i in mot_100.values()]
# p_m_100 = shapiro(x_m_100)
# print(f'Коэффициент Шапиро-Уилка для мотивов 100 ортологов: {p_m_100[0]}')
# print(f'Уровень значимости: {p_m_100[1]}')


# x_m_10 = [float(i) for i in mot_10.values()]
# p_m_10 = shapiro(x_m_10)
# print(f'Коэффициент Шапиро-Уилка для мотивов 10 ортологов: {p_m_10[0]}')
# print(f'Уровень значимости: {p_m_10[1]}')


# Построение графика 
nucl = ['a', 't', 'c', 'g']
acids = ['F', 'S', 'Y', 'C', 'L', 'W', 'P', 'H', 'R', 'Q', 'I', 'T', 'N',
         'K', 'M', 'V', 'A', 'D', 'G', 'E']
colors = ['#AE59DF', '#DF6B59', '#8ADF59', '#59CDDF']

nucleotides = your_gene.ntseq[0]
aa = your_gene.aaseq[0]

pie_nucl = [nucleotides.count(i) for i in nucl] 
bar_aa = [aa.count(i) for i in acids]

fig = plt.figure()
ax_1 = fig.add_subplot(2, 2, 1)
ax_2 = fig.add_subplot(2, 2, 2)

ax_1.pie(pie_nucl, labels=nucl, colors=colors, radius=1.2, autopct=f"{'%.1f%'}%", wedgeprops={'linewidth':1, 'edgecolor':'k'})
ax_2.bar(acids, bar_aa, 0.6, color='#AE59DF')
ax_1.set_title('Содержание нуклеотидов', fontsize=10)
ax_2.set_title('Содержание аминокислот')
ax_2.set_ylabel('Количество', fontsize=10)
plt.xticks(fontsize=7.5)
plt.yticks(fontsize=7.5)
plt.show()

