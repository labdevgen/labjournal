# Неизвестный палиндром

В результатах Exo-C с помощью приложения *FastQC* был обнаружен неизвестный палиндромный участок `CTCAGCGCTGAG`.
Частота его встречаемости, согласно данным скрипта *palindrome.py*, составляет 22.69%, в начале (5' конец, вторая позиция) - 6.88%.

```python
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
```

Было решено узнать, что собой представляет этот участок.

## Основная гипотеза

Палиндром является сдвоенным фрагментом blunt-адаптера, по какой-то причине потерявшим конец.

## Ход работы

С помощью *cutadapt* были найдены последовательности, содержащие искомый палиндром.

```
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
```

Далее результаты были снова обработаны с помощью *FastQC* ([данные здесь](http://htmlpreview.github.io/?https://github.com/regnveig/labjournal/blob/master/FastQC_results/fastqc_190718_1327.html)).
Программа обнаружила 3 длинных оверрепрезентированных последовательности, две из которых были определены как **TruSeq Adapter**, и одну неизвестную последовательность.
Их частоты составляют 0.2% (для TruSeq-адаптеров) и 0.1% (неизвестная последовательность) среди ридов, содержащих искомый палиндром.

Далее последовательности были выровнены относительно друг друга.

```
-----------------------------GAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTAT
------------------GGATCCCTCAGCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCA
GGATCCCTCAGCGCTGAGGGATCCCTCAGCAGATCGGAAGAGCACACGTC
```

Искомый палиндром был обнаружен в третьей (неизвестной) последовательности - GGATCC**CTCAGCGCTGAG**GGATCCCTCAGCAGATCGGAAGAGCACACGTC.
Также выяснено, что палиндром является частью ещё более крупного палиндрома, входящего в эту последовательность - GGATCCCTCAGCGCTGAGGGATCC.

В сочетании с вырезанным нами палиндромом он даёт ещё более длинный палиндром - **CCTCAGC***GCTGAGG*GATC**CCTCAGC***GCTGAGG*GATC**CCTCAGC**AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTAT

Было решено построить детально модель HiC и секвенирования по методике, чтобы понять, что произошло.

## Модель HiC и секвенирования

Методика описана в следующих статьях: [принципиальная схема](https://www.ncbi.nlm.nih.gov/pubmed/25437436), [подробная инструкция](https://www.ncbi.nlm.nih.gov/pubmed/29382556).

1. Фиксация, выделение ядер
2. Разрезание хроматина ДНКазой I
3. Обработка ДНК-полимеразой, фрагментом Кленова (5'→3' полимеразная активность, корректорная 3'→5' кусь-активность)
4. dA-tailing - фрагмент Кленова, dATP
5. ДНК-лигаза, Т-tailed биотин-меченый bridge-адаптер, blunt-ended Bridge безбиотиновый

```
bridge                      blunt

 P        Biot
 |        |
 5-GCTGAGGGATC-3           5-GCTGAGGGAC-3
3-TCGACTCC-5               3-CGACTCC-5

reversed 
         :egdirB             :tnulB
    5-CCTCAGCT-3             5-CCTCAGC-3
3-CTAGGGAGTCG-5           3-CAGGGAGTCG-5
     |      |
     Biot   P
     
     
Products of adapter ligation:

(genome-A-)Bridge/GATC/egdirB-(T-genome):
        (A)GCTGAGG/GATC/CCTCAGC(T)
just sequence, palyndromic: AGCTGAGGGATCCCTCAGCT

(genome-A-)Bridge---egdirB(T-genome):
        (A)GCTGAGG-CCTCAGC(T)
just sequence, palyndromic: AGCTGAGGCCTCAGCT

(remove GATC)egdirB(remove T)-(remove T)Bridge(remove GATC):
             CCTCAGC-------------------GCTGAGG
just sequence, palyndromic: CCTCAGCGCTGAGG

       tnulB-Blunt:
(GTC)CCTCAGC-GCTGAGG(GAC)
Just sequence (palyndromic): GTCCCTCAGCGCTGAGGGAC

```


6. Полинуклеотидкиназа (прикрепляет фосфат к 5'), затем лигаза
7. Растворение белков и очистка ДНК
8. ДНК-полимераза, dATP, dGTP - достраивание цепей
9. Фрагментация ДНК до размеров 100-300

Далее методика меняется на протокол NEBNext Ultra II ([ссылка](http://www.bea.ki.se/documents/datasheet_NEB_Ultra%20II%20DNA.pdf)).

...

**Гипотезы**
1. Наш палиндром `CTCAGCGCTGAG` (23% ридов) может быть либо стыком 
egdirB(remove T)-(remove T)Bridge

либо

tnulB-Blunt

Различить эти два сценария можно только анализируя буквы до/после найденного палиндрома. А в случае, если они "обкусаны", вообще нельзя.

2. Длинная оверепрезентированная последовательность (**CCTCAGC***GCTGAGG*GATC**CCTCAGC***GCTGAGG*GATC**CCTCAGC**AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTAT) набиолее вероятно обарзовалась так: сначала слиплись две последовательности Bridge, образовав egdirB(remove T)-(remove T)Bridge-GATC. Затем, они слиплись друг с другом по липким GATC концам.

Т.е. сначала бриджи слиплись друг с другом "спинками", а потом два таких соеденились по липкому концу.

egdirB(remove T)-(remove T)Bridge-GATC-egdirB(remove T)-(remove T)Bridge

## Промежуточная задача

Поискать bridge и blunt-адаптеры в библиотеке.
Bridge - `AGCTGAGGGATC`, blunt - `GCTGAGGGAC` и палиндром-содержащий участок `CCTCAGCGCTGAGGGAC`.
Найти их общее количество, а также распределение по позициям в риде.

Обработка была произведена с помощью скрипта `palindrome2.py`.

```python
import gzip
import sys
import pandas as pd
import string

def Out(found_, total_):
    print("Found: %d | Total: %d (%4f%%)" % (found_, total_, found_ * 100 / total_), end='\r')

filename = './sample.fastq.gz'
seq = 'CCTCAGCGCTGAGGGAC'
comment = "Double Blunt adapter"

print(f"\nHi there.\nWe're looking for: {seq} ({comment})\n")

input0 = gzip.open(filename, 'r')

counter = 0
total = 1
found = 0
tyk = 0
df = pd.DataFrame({
'position' : [0],
'count' : [0]
})

for line in input0:

    Out(found, total)

    counter += 1
    if counter == 5:
        counter = 1
    if counter < 1:
        continue

    if (counter == 2):
        tyk = line.decode("utf-8").find(seq)
        if (tyk == -1):
            counter = -2
            total += 1
            continue


    if counter == 4:
        if df.loc[df['position'] == tyk].empty:
            df = df.append({'position': tyk, 'count': 1}, ignore_index=True)
        else:
            df.at[df.loc[df['position'] == tyk].index[0], 'count'] += 1

        found += 1
        total += 1

output0 = open('./report_' + seq + '.txt', 'w')
df.sort_values(by=['position'], ascending=True).to_string(output0)
print('\n')

input0.close()
output0.close()
```

Результаты:

* Двойной blunt-адаптер встречается в 0.6% ридов.
Пики наблюдаются на позициях 0, 11, небольшой пик на 18 ([данные](./scripts_results/report_palindrome_doubleblunt_190719.txt)).

```
$ python3 ./palindrome_doubleblunt.py

Hi there.
We're looking for: CCTCAGCGCTGAGGGAC (Double Blunt adapter)

Found: 775465 | Total: 128195238 (0.604909%)
```

![График doubleblunt](./scripts_results/graph_doubleblunt_190719.png)

* Одиночный blunt-адаптер встречается в 1.42% ридов.
Пики на позициях 0, 7 и 18 ([данные](./scripts_results/report_palindrome_blunt_190719.txt)).

```
$ python3 ./palindrome_blunt.py

Hi there.
We're looking for: GCTGAGGGAC (Blunt adapter)

Found: 1816057 | Total: 128195238 (1.416634%)
```

![График blunt](./scripts_results/graph_blunt_190719.png)

**TODO:**

поискать ешё димеры egdirb-bridge-GATC

Или, что может быть ещё лучше, проанализировать контент букв до и после *CTCAGCGCTGAG*
