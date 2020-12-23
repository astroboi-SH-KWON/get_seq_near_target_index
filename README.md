# get_seq_near_target_index


[Main_SY.py 관련](./Main_SY.py)

sgRNA design for selective knockout of dominant mutant alleles
(guide sequence에는 mismatch가 있어도 activity가 어느 정도 유지될 수 있으나, PAM에 위치할 경우 activity가 거의 나타나지 않음을 이용)

1. ClinVar data에서 dominant mutation만 filter함


2. Filtered ClinVar data에서 mutation에 의해 새로 생긴 PAM을 기준으로 sgRNA와 주변 context를 가져옴
    2.1. Reference sequence: filtered ClinVar data의 1st (#CHROM), 2nd (POS) colum을 참고하여 reference human genome (GRCh38_p13)에서 position을 기준으로 해당 gene의 sequence를 가져옴 (기존의 mutation position 기준 양쪽 60 bp를 가장 큰 gene의 기준인 2.3 Mbp로 바꿀 수 있으나.. 너무 비효율적인 방법인 듯하여 이렇게 정하였어요..)
    2.2. Mutation sequence: reference sequence에서 filter ClinVar data의 3rd column (REF)의 sequence를 4th column (ALT)의 sequence로 치환
    2.3. Mutation sequence 상에서 ALT sequence에 overlapping되는 PAM sequence를 전부 찾고 주변의 context를 가져옴 (중복은 아래에서 한번에 제거됨)
    2.4. 2.3.의 context에서 guide를 가져와 만약 2.1.의 해당 gene 내에서 가능한 guide (PAM의 존재에 따라)의 목록과 비교하여 해당하면 제거

fig.1
	         <u>|5’ context	|Guide	|PAM	    |3’ context</u>	       
SaCas9	        |22 nt	    |21 nt	|NNGRRT	    |3 nt
SaCas9-KKH		|           |21 nt	|NNNRRT	    |
SaCas9-NNG		|           |21 nt	|NNG	    |
SauriCas9		|           |21 nt	|NNGG	    |
SauriCas9-KKH	|           |21 nt	|NNRG	    |
St1Cas9		    |           |19 nt	|NNRGAA	    |
Nm1Cas9		    |           |23 nt	|NNNNGATT	|
Nm2Cas9		    |           |22 nt	|NNNNCC	    |
CjCas9		    |           |22 nt	|NNNNRYAC   |

fig.2 Example, SaCas9-NNG (PAM: 5’-NNG-3’)
![alt text](./fig_2.PNG)

Arrow는 guide + PAM을 나타냄

Blue, orange, red (같은 gene 내의 다른 부위에 PAM을 가지는 같은 guide 서열이 존재)는 탈락
따라서, black만 통과




+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
input : 
    ALL_CDS_INFO = "all_ccds_filtered_201130_CCDS_" + TYPE + "_current.txt"  # 20201217
    MUT_INFO = [ClinVar_data : ./input/200907_Dominant filter.txt](./input/200907_Dominant filter.txt)
    
 
1. get rows in ClinVar_data(= MUT_INFO) if POS is in CDS(ALL_CDS_INFO) 
    1-1. get seq pairs with WIN_SIZE near target seq(= POS column in ClinVar_data, len(REF column)) in CDS not in whole genome
    



??????????    
len(ref in cds) = 100
60 bp + ALT + 60 bp > 120

len(ref in cds) = 500
60 bp + ALT + 60 bp > 120
=> should start after 121 in (ref in cds)




