# Generating Reed-Solomon codebooks

This folder contains multiple notebooks for generating Reed-Solomon codebooks for
use with seqFISH experiments. Different notebooks construct different kinds of Reed-Solomon codes. Some generate Reed-Solomon codes over a base field (prime number of elements), and others generate Reed-Solomon codes over an extension field (prime power of elements with integer power > 1), extended Reed-Solomon codes (codes with n = q or n = q+1).


To see what characteristics are possible to get from different Reed-Solomon codebooks, you can first look at [Make_RS_Code_Table.jl.ipynb](https://github.com/CaiGroup/UntanglingBarcodes/blob/main/codebook_generation/get_RS_codebooks/Make_RS_Code_Table_extended.jl.ipynb) ([Google Colab](https://colab.research.google.com/github/CaiGroup/UntanglingBarcodes/blob/master/codebook_generation/get_RS_codebooks/colab/Make_RS_Code_Table.jl.ipynb) which produces the table of Reed-Solomon codes shown at the bottom of this page. Codes are sorted by their relevant characteristics including n (the number of barcoding blocks, q (the number of pseudocolors), minimum Hamming distance (between any two codewords), codeword weight enumeration (the number of codewords in the code that with a given number of non-zero symbols). The number of images in a seqFISH implementation of a code is n*q. This table will give help gives a guide to choosing a code with appropriate characteracterists for designing a given experiment.



* [gen_RS_q11k7_code.jl.ipynb](https://github.com/CaiGroup/UntanglingBarcodes/blob/main/codebook_generation/get_RS_codebooks/gen_RS_q11k7_code.jl.ipynb) ([Google Colab](https://colab.research.google.com) generates the Reed-Solomon code implemented in experiment. See processing workflow [here](https://github.com/CaiGroup/UntanglingBarcodes/tree/main/real_data_processing_workflows/Reed-Solomon_encoded_experiment_workflow). The notebook can also generate tables of Reed-Solomon codewords over base finite fields (where q is prime).
 * [gen_extended_codes_cuda.jl.ipynb](https://github.com/CaiGroup/UntanglingBarcodes/blob/main/codebook_generation/get_RS_codebooks/gen_extended_codes_cuda.jl.ipynb) finds codewords in extended Reed-Solomon codes which have n = 1 or n = q+1. This notebook uses cuda GPUs since as many matrix multiplications are required.
 * [get_shortened_cb.jl.ipynb](https://github.com/CaiGroup/UntanglingBarcodes/blob/main/codebook_generation/get_RS_codebooks/get_shortened_cb.jl.ipynb) Finds shortened Reed-Solomon codes with n < q-1.
 * [gen_extended_codes_ext_fields_q9_n9-10_nmk4.ipynb](https://github.com/CaiGroup/UntanglingBarcodes/blob/main/codebook_generation/get_RS_codebooks/gen_extended_codes_ext_fields_q9_n9-10_nmk4.ipynb) ([Google Colab](https://colab.research.google.com/github/CaiGroup/UntanglingBarcodes/blob/master/codebook_generation/get_RS_codebooks/colab/gen_extended_codes_ext_fields_q9_n9-10_nmk4.ipynb)) generates Reed-Solomon codes over extended finite fields (where q is a prime power greater than 1). This notebook also generates extended codes (where n = q or n=q+1) over the same finite fields. I do not provide GPU acceleration for the extended fields here since the most practical Reed-Solomon codes over extended fields for seqFISH have small qs of 8 or 9 allowing them to be computed quickly on a CPU.


The Manifest.toml and Project.toml files contain the Julia environment to run these notebooks.


|# of Images|q  |n  |mhd|w3    |w4      |w5       |w6         |
|-----|---|---|---|------|--------|---------|-----------|
|16   |5  |4  |2  |48    |52      |0        |0          |
|20   |5  |5  |2  |120   |260     |204      |0          |
|24   |5  |6  |2  |240   |780     |1224     |820        |
|20   |5  |5  |3  |40    |40      |44       |0          |
|24   |5  |6  |3  |80    |120     |264      |160        |
|24   |5  |6  |4  |0     |60      |24       |40         |
|24   |7  |4  |2  |120   |186     |0        |0          |
|30   |7  |5  |2  |300   |930     |1110     |0          |
|36   |7  |6  |2  |600   |2790    |6660     |6666       |
|42   |7  |7  |2  |1050  |6510    |23310    |46662      |
|48   |7  |8  |2  |1680  |13020   |62160    |186648     |
|24   |7  |4  |3  |24    |24      |0        |0          |
|30   |7  |5  |3  |60    |120     |162      |0          |
|36   |7  |6  |3  |120   |360     |972      |948        |
|42   |7  |7  |3  |210   |840     |3402     |6636       |
|48   |7  |8  |3  |336   |1680    |9072     |26544      |
|36   |7  |6  |4  |0     |90      |108      |144        |
|42   |7  |7  |4  |0     |210     |378      |1008       |
|48   |7  |8  |4  |0     |420     |1008     |4032       |
|42   |7  |7  |5  |0     |0       |126      |84         |
|48   |7  |8  |5  |0     |0       |336      |336        |
|28   |8  |4  |2  |168   |301     |0        |0          |
|35   |8  |5  |2  |420   |1505    |2100     |0          |
|42   |8  |6  |2  |840   |4515    |12600    |14707      |
|49   |8  |7  |2  |1470  |10535   |44100    |102949     |
|56   |8  |8  |2  |2352  |21070   |117600   |411796     |
|63   |8  |9  |2  |3528  |37926   |264600   |1235388    |
|35   |8  |5  |3  |70    |175     |266      |0          |
|42   |8  |6  |3  |140   |525     |1596     |1834       |
|49   |8  |7  |3  |245   |1225    |5586     |12838      |
|56   |8  |8  |3  |392   |2450    |14896    |51352      |
|63   |8  |9  |3  |588   |4410    |33516    |154056     |
|42   |8  |6  |4  |0     |105     |168      |238        |
|49   |8  |7  |4  |0     |245     |588      |1666       |
|56   |8  |8  |4  |0     |490     |1568     |6664       |
|63   |8  |9  |4  |0     |882     |3528     |19992      |
|49   |8  |7  |5  |0     |0       |147      |147        |
|56   |8  |8  |5  |0     |0       |392      |588        |
|63   |8  |9  |5  |0     |0       |882      |1764       |
|40   |9  |5  |2  |560   |2280    |3640     |0          |
|48   |9  |6  |2  |1120  |6840    |21840    |29128      |
|56   |9  |7  |2  |1960  |15960   |76440    |203896     |
|64   |9  |8  |2  |3136  |31920   |203840   |815584     |
|72   |9  |9  |2  |4704  |57456   |458640   |2446752    |
|80   |9  |10 |2  |6720  |95760   |917280   |6116880    |
|40   |9  |5  |3  |80    |240     |408      |0          |
|48   |9  |6  |3  |160   |720     |2448     |3232       |
|56   |9  |7  |3  |280   |1680    |8568     |22624      |
|64   |9  |8  |3  |448   |3360    |22848    |90496      |
|72   |9  |9  |3  |672   |6048    |51408    |271488     |
|80   |9  |10 |3  |960   |10080   |102816   |678720     |
|48   |9  |6  |4  |0     |120     |240      |368        |
|56   |9  |7  |4  |0     |280     |840      |2576       |
|64   |9  |8  |4  |0     |560     |2240     |10304      |
|72   |9  |9  |4  |0     |1008    |5040     |30912      |
|80   |9  |10 |4  |0     |1680    |10080    |77280      |
|56   |9  |7  |5  |0     |0       |168      |224        |
|64   |9  |8  |5  |0     |0       |448      |896        |
|72   |9  |9  |5  |0     |0       |1008     |2688       |
|80   |9  |10 |5  |0     |0       |2016     |6720       |
|70   |11 |7  |2  |3150  |31850   |190890   |636370     |
|80   |11 |8  |2  |5040  |63700   |509040   |2545480    |
|90   |11 |9  |2  |7560  |114660  |1145340  |7636440    |
|100  |11 |10 |2  |10800 |191100  |2290680  |19091100   |
|110  |11 |11 |2  |14850 |300300  |4199580  |42000420   |
|120  |11 |12 |2  |19800 |450450  |7199280  |84000840   |
|70   |11 |7  |3  |350   |2800    |17430    |57820      |
|80   |11 |8  |3  |560   |5600    |46480    |231280     |
|90   |11 |9  |3  |840   |10080   |104580   |693840     |
|100  |11 |10 |3  |1200  |16800   |209160   |1734600    |
|110  |11 |11 |3  |1650  |26400   |383460   |3816120    |
|120  |11 |12 |3  |2200  |39600   |657360   |7632240    |
|70   |11 |7  |4  |0     |350     |1470     |5320       |
|80   |11 |8  |4  |0     |700     |3920     |21280      |
|90   |11 |9  |4  |0     |1260    |8820     |63840      |
|100  |11 |10 |4  |0     |2100    |17640    |159600     |
|110  |11 |11 |4  |0     |3300    |32340    |351120     |
|120  |11 |12 |4  |0     |4950    |55440    |702240     |
|70   |11 |7  |5  |0     |0       |210      |420        |
|80   |11 |8  |5  |0     |0       |560      |1680       |
|90   |11 |9  |5  |0     |0       |1260     |5040       |
|100  |11 |10 |5  |0     |0       |2520     |12600      |
|110  |11 |11 |5  |0     |0       |4620     |27720      |
|120  |11 |12 |5  |0     |0       |7920     |55440      |
|108  |13 |9  |2  |11088 |201096  |2411640  |19294128   |
|120  |13 |10 |2  |15840 |335160  |4823280  |48235320   |
|132  |13 |11 |2  |21780 |526680  |8842680  |106117704  |
|144  |13 |12 |2  |29040 |790020  |15158880 |212235408  |
|156  |13 |13 |2  |37752 |1141140 |24633180 |394151472  |
|168  |13 |14 |2  |48048 |1597596 |38318280 |689765076  |
|108  |13 |9  |3  |1008  |15120   |185976   |1483776    |
|120  |13 |10 |3  |1440  |25200   |371952   |3709440    |
|132  |13 |11 |3  |1980  |39600   |681912   |8160768    |
|144  |13 |12 |3  |2640  |59400   |1168992  |16321536   |
|156  |13 |13 |3  |3432  |85800   |1899612  |30311424   |
|168  |13 |14 |3  |4368  |120120  |2954952  |53044992   |
|108  |13 |9  |4  |0     |1512    |13608    |114912     |
|120  |13 |10 |4  |0     |2520    |27216    |287280     |
|132  |13 |11 |4  |0     |3960    |49896    |632016     |
|144  |13 |12 |4  |0     |5940    |85536    |1264032    |
|156  |13 |13 |4  |0     |8580    |138996   |2347488    |
|168  |13 |14 |4  |0     |12012   |216216   |4108104    |
|108  |13 |9  |5  |0     |0       |1512     |8064       |
|120  |13 |10 |5  |0     |0       |3024     |20160      |
|132  |13 |11 |5  |0     |0       |5544     |44352      |
|144  |13 |12 |5  |0     |0       |9504     |88704      |
|156  |13 |13 |5  |0     |0       |15444    |164736     |
|168  |13 |14 |5  |0     |0       |24024    |288288     |
|180  |16 |12 |2  |46200 |1566675 |37588320 |657809460  |
|195  |16 |13 |2  |60060 |2262975 |61081020 |1221646140 |
|210  |16 |14 |2  |76440 |3168165 |95014920 |2137880745 |
|225  |16 |15 |2  |95550 |4320225 |142522380|3563134575 |
|240  |16 |16 |2  |117600|5760300 |207305280|5701015320 |
|255  |16 |17 |2  |142800|7532700 |293682480|8810660040 |
|180  |16 |12 |3  |3300  |96525   |2352240  |41108760   |
|195  |16 |13 |3  |4290  |139425  |3822390  |76344840   |
|210  |16 |14 |3  |5460  |195195  |5945940  |133603470  |
|225  |16 |15 |3  |6825  |266175  |8918910  |222672450  |
|240  |16 |16 |3  |8400  |354900  |12972960 |356275920  |
|255  |16 |17 |3  |10200 |464100  |18378360 |550608240  |
|180  |16 |12 |4  |0     |7425    |142560   |2577960    |
|195  |16 |13 |4  |0     |10725   |231660   |4787640    |
|210  |16 |14 |4  |0     |15015   |360360   |8378370    |
|225  |16 |15 |4  |0     |20475   |540540   |13963950   |
|240  |16 |16 |4  |0     |27300   |786240   |22342320   |
|255  |16 |17 |4  |0     |35700   |1113840  |34529040   |
|180  |16 |12 |5  |0     |0       |11880    |152460     |
|195  |16 |13 |5  |0     |0       |19305    |283140     |
|210  |16 |14 |5  |0     |0       |30030    |495495     |
|225  |16 |15 |5  |0     |0       |45045    |825825     |
|240  |16 |16 |5  |0     |0       |65520    |1321320    |
|255  |16 |17 |5  |0     |0       |92820    |2042040    |
|208  |17 |13 |2  |68640 |2757040 |79382160 |1693513536 |
|224  |17 |14 |2  |87360 |3859856 |123483360|2963648688 |
|240  |17 |15 |2  |109200|5263440 |185225040|4939414480 |
|256  |17 |16 |2  |134400|7017920 |269418240|7903063168 |
|272  |17 |17 |2  |163200|9177280 |381675840|12213824896|
|288  |17 |18 |2  |195840|11799360|528474240|18320737344|
|208  |17 |13 |3  |4576  |160160  |4674384  |99610368   |
|224  |17 |14 |3  |5824  |224224  |7271264  |174318144  |
|240  |17 |15 |3  |7280  |305760  |10906896 |290530240  |
|256  |17 |16 |3  |8960  |407680  |15864576 |464848384  |
|272  |17 |17 |3  |10880 |533120  |22474816 |718402048  |
|288  |17 |18 |3  |13056 |685440  |31118976 |1077603072 |
|208  |17 |13 |4  |0     |11440   |267696   |5875584    |
|224  |17 |14 |4  |0     |16016   |416416   |10282272   |
|240  |17 |15 |4  |0     |21840   |624624   |17137120   |
|256  |17 |16 |4  |0     |29120   |908544   |27419392   |
|272  |17 |17 |4  |0     |38080   |1287104  |42375424   |
|288  |17 |18 |4  |0     |48960   |1782144  |63563136   |
|208  |17 |13 |5  |0     |0       |20592    |329472     |
|224  |17 |14 |5  |0     |0       |32032    |576576     |
|240  |17 |15 |5  |0     |0       |48048    |960960     |
|256  |17 |16 |5  |0     |0       |69888    |1537536    |
|272  |17 |17 |5  |0     |0       |99008    |2376192    |
|288  |17 |18 |5  |0     |0       |137088   |3564288    |
