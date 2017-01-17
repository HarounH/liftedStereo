all:
	g++ hoStereo4.cpp -std=c++11 -lpng -I ../libs/opengm-master/src/external/QPBO-v1.3.src-patched/ -lexternal-library-qpbo-shared -o a.out -Wfatal-errors
a:
	g++ stereoInfer.cpp -std=c++11 -lpng -Wfatal-errors -o a.out
b:
	g++ stereoInfer2.cpp -std=c++11 -lpng -Wfatal-errors -o b.out
c:
	g++ stereoInfer3.cpp -std=c++11 -lpng -Wfatal-errors -o c.out	
testc:
	g++ helloOpenGM.cpp  -std=c++11 -lpng -Wfatal-errors -o t.out
biga:
	./a.out single ../../data/tsubaka/scene1.row3.col3.png ../../data/tsubaka/scene1.row3.col5.png ../../output/tsu_biga.png 30 50 3 200 ../../data/tsubaka/truedisp.row3.col3.png
bigb:
	./b.out single ../../data/tsubaka/scene1.row3.col3.png ../../data/tsubaka/scene1.row3.col5.png ../../output/tsu_bigb.png 30 50 3 200 2 10.0 ../../data/tsubaka/truedisp.row3.col3.png ../../data/tsubaka/symsl1.png
bigc:
	./c.out single ../../data/tsubaka/scene1.row3.col3.png ../../data/tsubaka/scene1.row3.col5.png ../../output/tsu_bigc.png 30 50 3 200 2 64 ../../data/tsubaka/truedisp.row3.col3.png
bigall: a b c
	./a.out single ../../data/tsubaka/scene1.row3.col3.png ../../data/tsubaka/scene1.row3.col5.png ../../output/tsu_biga.png 30 50 3 200 ../../data/tsubaka/truedisp.row3.col3.png > ../../output/biga_20_100_3_200.stats
	./b.out single ../../data/tsubaka/scene1.row3.col3.png ../../data/tsubaka/scene1.row3.col5.png ../../output/tsu_bigb.png 30 50 3 200 2 30.0 ../../data/tsubaka/truedisp.row3.col3.png > ../../output/bigb_20_100_3_200_2_1000.stats
	./c.out single ../../data/tsubaka/scene1.row3.col3.png ../../data/tsubaka/scene1.row3.col5.png ../../output/tsu_bigc.png 30 50 3 200 2 64 ../../data/tsubaka/truedisp.row3.col3.png	> ../../output/bigb_20_100_3_200_2_64.stats

newa:
	./a.out single ../../data/all2/Art/view1.png ../../data/all2/Art/view5.png ../../data/all2/Art/my_disp1.png 36 25 3 200 ../../data/all2/Art/disp1.png
newb:
	./b.out single ../../data/all2/Art/view1.png ../../data/all2/Art/view5.png ../../data/all2/Art/l1_disp.png 30 50 3 200 2 10.0 ../../data/all2/Art/disp1.png ../../data/all2/Art/symsl1.png
testr:
	./t.out
clean:
	rm *.out