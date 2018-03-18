# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 10:57:49 2018

@author: Darida
"""

from Bio import SeqIO
#Класс, описывающий к-меры. Имеет характеристики:счётчик,позиция,№хромосомы, для которых определены соответствующие функции: увеличение счётчика, присвоение
class Kmer:
    sequence = ''    

    def __init__(self,kmer_name):
        self.sequence = kmer_name
        self.counter = 1
        self.lst = []
        self.chr = []
    
    def increase(self):
        self.counter += 1
    
    def wr_N_chr(self, i):
        self.chr.append(i)
    def add_locus(self, i):
        self.lst.append(i)
#Адреса входного и выходного файлов 
in_file = 'C:/Users/Darida/Documents/seq_y_pestis.fasta'
out_file =  'C:/Users/Darida/Documents/seq_y_pestis_out.txt'
#Чтение последовательности из файла
with open(in_file, "rU") as handle:
    seq=[]   
    for record in SeqIO.parse(handle, "fasta"):
        seq.append(str(record.seq))

seq_n = len(seq)
seq_lng = len(record.seq)
print(seq_lng)
print(seq_n)
#Длина к-мера
kmer_size = 23
kmer_dict = {}
#Для каждой хромосомы
for index_chr in range(seq_n):   
    #Проверяй к-меры на каждой позиции
    for index in range(seq_lng-kmer_size+1): 
        current_kmer = seq[index_chr][index:(index+kmer_size)] 
        #Если такой уже есть в словаре, увеличь счётчик и добавь новый локус
        if current_kmer in kmer_dict:
            kmer_dict[current_kmer].increase()
            kmer_dict[current_kmer].add_locus(index)
        #Если нет, добавь новый, запиши его № хромосомы и текущий локус
        else:
            kmer_dict[current_kmer] = Kmer(current_kmer)
            kmer_dict[current_kmer].add_locus(index)
            kmer_dict[current_kmer].wr_N_chr(index_chr)
             
max = 0
#Запись в выходной файл
g = open(out_file, "w")
#Строка с названиями колонок
g.write("k-мер"+"\t"+"\t"+"\t"+"Счётчик "+"\t"+"N хромосомы" +"\t"+"\t"+"Координаты"+"\n")

for k,s in kmer_dict.items():
    #Определения наиболее частого к-мера
    if s.counter > max:
        max = s.counter
    #Запись текущего словаря(к-мер, счётчки, № хромосомы, координаты) . закомментирован для ускорения
    #g.write(k+"\t"+"\t"+str(s.counter)+"\t"+str(s.chr)+"\t"+"\t"+"\t"+str(s.lst)+"\n")
#Вывод информации о наиболее частом к-мере
g.write("Местоположения наиболее частого k-мера\n")
for k,s in kmer_dict.items():
    if s.counter == max:
        g.write(k+"\t"+"\t"+str(s.counter)+"\t"+str(s.chr)+"\t"+"\t"+"\t"+str(s.lst)+"\n")
g.write("\n")
g.close()