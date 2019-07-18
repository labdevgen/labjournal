# Лабораторный журнал

*Автор:* Валеев Э.С.

## Неизвестный палиндром

В результатах Exo-C с помощью приложения *FastQC* был обнаружен неизвестный палиндромный участок `CTCAGCGCTGAG`.
Частота его встречаемости, согласно данным скрипта *palindrome.py*, составляет 22.69%, в начале (5' конец, вторая позиция) - 6.88%.

<details>
<summary>palindrome.py</summary>
<p><pre>
	import gzip
	import sys
	
	def Out(found_, total_):
	    print("Found: %d | Total: %d" % (found_, total_), end='\r')
	
	filename = './sample.fastq.gz'
	
	input0 = gzip.open(filename, 'r')
	output0 = open('./output.fq', 'w')
	
	new_found_land = ""
	counter = 0
	total = 0
	found = 0
	
	for line in input0:
	
	    Out(found, total)
	
	    counter += 1
	    if counter == 5:
	        counter = 1
	    if counter < 1:
	        continue
	
	    tyk = line.decode().find("CTCAGCGCTGAG")
	
	    if (counter == 2) and (tyk == -1):
	        new_found_land = ""
	        counter = -2
	        total += 1
	        continue
	
	    new_found_land += line.decode("utf-8")
	
	    if counter == 4:
	        output0.write(new_found_land)
	        new_found_land = ""
	        found += 1
	        total += 1
	
	print('\n')
	
	input0.close()
	output0.close()
</pre></p>
</details>

> Hello
>> Hellon't

Было решено узнать, что собой представляет этот участок.

### Основная гипотеза

Палиндром является сдвоенным фрагментом blunt-адаптера `GCTGAGGGAC`, по какой-то причине потерявшим `GGAC`.

### Ход работы

С помощью *cutadapt* были найдены последовательности, содержащие искомый палиндром.

	$ cutadapt -g ^CCTCAGCGCTGAG --trimmed-only -o ./output.ca.fastq ./sample.fastq.gz
	This is cutadapt 1.18 with Python 3.7.3
	Command line parameters: -g ^CCTCAGCGCTGAG --trimmed-only -o ./output.ca.fastq ./sample.fastq.gz
	Processing reads on 1 core in single-end mode ...
	Finished in 1481.04 s (12 us/read; 5.19 M reads/minute).
	
	=== Summary ===
	
	Total reads processed:             128,195,237
	Reads with adapters:                10,090,133 (7.9%)
	Reads written (passing filters):    10,090,133 (7.9%)
	
	Total basepairs processed: 19,229,285,550 bp
	Total written (filtered):  1,383,201,351 bp (7.2%)
	
	=== Adapter 1 ===
	
	Sequence: CCTCAGCGCTGAG; Type: anchored 5'; Length: 13; Trimmed: 10090133 times.
	
	No. of allowed errors:
	0-9 bp: 0; 10-13 bp: 1
	
	Overview of removed sequences
	length  count   expect  max.err error counts
	12      938313  7.6     1       0 938313
	13      9066637 1.9     1       8892621 174016
	14      85183   1.9     1       0 85183

Далее результаты были снова обработаны с помощью *FastQC* ([данные здесь](./FastQC_results/fastqc_190718_1327.html)).
Программа обнаружила 3 длинных оверрепрезентированных последовательности, две из которых были определены как **TruSeq Adapter**, и одну неизвестную последовательность.
Их частоты составляют 0.2% (для TruSeq-адаптеров) и 0.1% (неизвестная последовательность) среди ридов, содержащих искомый палиндром.

Далее последовательности были выровнены относительно друг друга.

	-----------------------------GAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTAT
	------------------GGATCCCTCAGCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCA
	GGATCCCTCAGCGCTGAGGGATCCCTCAGCAGATCGGAAGAGCACACGTC

Искомый палиндром был обнаружен в третьей (неизвестной) последовательности - GGATCC**CTCAGCGCTGAG**GGATCCCTCAGCAGATCGGAAGAGCACACGTC.
Также выяснено, что палиндром является частью ещё более крупного палиндрома, входящего в эту последовательность - GGATCCCTCAGCGCTGAGGGATCC.
В сочетании с вырезанным нами палиндромом он даёт ещё более длинный палиндром - CCTCAGC**GCTGAGGGATCCCTCAGCGCTGAGGGATCCCTCAGC**AGATCGGAAGAGCACACGTC.
Было решено построить детально модель HiC и секвенирования по методике, чтобы понять, что произошло.

### Модель HiC и секвенирования

Методика и принцип описаны в следующих статьях: [раз](https://www.ncbi.nlm.nih.gov/pubmed/25437436), [два](https://www.ncbi.nlm.nih.gov/pubmed/29382556).

1. Фиксация, выделение ядер
2. Разрезание хроматина ДНКазой I
3. Обработка ДНК-полимеразой, фрагментом Кленова (5'-3' полимеразная активность, корректорная 3'-5' кусь-активность)
4. dA-tailing - фрагмент Кленова, dATP
5. ДНК-лигаза, Т-tailed биотин-меченый bridge-адаптер, blunt-ended Bridge безбиотиновый

**Структуры:**

	bridge                 blunt
	
	 P        Biot
	 |        |
	 GCTGAGGGATC           GCTGAGGGAC
	TCGACTCC               CGACTCC

6. Полинуклеотидкиназа (прикрепляет фосфат к 5'), затем лигаза
7. Растворение белков и очистка ДНК
8. ДНК-полимераза, dATP, dGTP
