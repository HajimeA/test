#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//一つの世代を構成する個体数 (交叉の条件上、偶数にすること)
#define IND_NUM 20
//遺伝子の長さ　長すぎると計算中にオーバーフローするので注意
#define GENE_LEN 16
//N進法（ここには指定したい進法から1つ引いた数にすること　大きすぎると計算中にオーバーフローするので注意）
#define ARY 1
//世代交代の回数
#define SEDAI_MAX 50
//交叉率
#define CROS_RATE 0.95
//突然変異率
#define MUTA_RATE 0.05

int Random(int min, int max)			//乱数発生関数 整数Ｎを受け取り０から(Ｎ－１)までの整数をランダムに返す
{
	return min + (int)(rand()*(max - min + 1.0) / (1.0 + RAND_MAX));
}

int N_ary_to_decimal(int gene_sanmple[GENE_LEN])		//N進法から10進法に変換する関数
{
	int decimal = 0;	//10進法
	int digit;			//桁数
	int base = 1;		//底

	for (digit = GENE_LEN; digit >0; digit--) {

		decimal = decimal + gene_sanmple[digit - 1] * base;
		base = base * (ARY + 1);
	}
	return decimal;
}

void main() {
	int i, j;						//計算用変数
	int counter;					//世代交代数

	int gene[IND_NUM][GENE_LEN];	//個体群の配列
	int gene_sanple[GENE_LEN];		//1つ個体の遺伝子列（計算用）
	int Phenotype[IND_NUM];			//個体の数値

	int elite_Ind;					//エリート個体の番号
	int elite_num;					//エリート個体の数値
	int elite[GENE_LEN];			//エリート個体の保存用配列

	double total, rulet_value, rulet_allow;	//ルーレット選択用変数
	int rulet_num;							//ルーレットで選択された個体の番号
	int temp_gene[IND_NUM][GENE_LEN];		//ルーレット選択で用いる個体群の退避用配列

	int cut;						//交叉における切断点を表す変数
	int temp;						//交叉における入替用変数
	double d_ret;					//交叉及び突然変異で用いる変数

	time_t t;						//乱数発生用の変数

	srand((unsigned)time(&t));	//time関数を利用して乱数の生成

	//個体群の生成
	for (i = 0; i < IND_NUM; i++) {			//個体数
		for (j = 0; j < GENE_LEN; j++) {	//遺伝子の長さ
			gene[i][j] = Random(0, ARY);	//各遺伝子にランダムに遺伝子を入力
		}

	}
	
	//生成された個体群の確認
	printf("生成された個体群\n");
	for (i = 0; i < IND_NUM; i++) {			//個体数
		printf("%03d:", i);
		for (j = 0; j < GENE_LEN; j++) {	//遺伝子の長さ
			printf("%d", gene[i][j]);		//各遺伝子の出力
			gene_sanple[j] = gene[i][j];
		}
		printf(":%d\n", N_ary_to_decimal(gene_sanple));
	}
	printf("\n");
	
	counter = 0;	//世代交代数をリセット

	while (counter < SEDAI_MAX) {		//世代交代数に達するまでループ

		//個体の評価　今回は10進数にするとより大きい程良いという条件とする。
		for (i = 0; i < IND_NUM; i++) {			//個体数分ループ
			for (j = 0; j < GENE_LEN; j++) {	//遺伝子の長さ分ループ
				gene_sanple[j] = gene[i][j];	//個体群の1つの遺伝子を配列に格納
			}	
			
			Phenotype[i] = N_ary_to_decimal(gene_sanple);	//各個体の数値を計算
		}

		//エリート個体の走査
		elite_Ind = 0;		//最初がエリートと仮定
		elite_num = Phenotype[0];			//一番最初の個体の成績を代入
		for (i = 1; i < IND_NUM; i++) {		//全ての個体の数値とエリートの数値を比べる
			if (elite_num < Phenotype[i]) {	//仮のエリートよりも成績のよいのが見つかったらエリート個体の番号と数値を更新する
				elite_Ind = i;
				elite_num = Phenotype[i];
			}
		}
		for (i = 0; i < GENE_LEN; i++) {	//elite_Indの遺伝子をエリート配列に入れる
			elite[i] = gene[elite_Ind][i];
		}


		//ルーレット選択
		total = 0.0;					//総和のリセット
		for (i = 0; i < IND_NUM; i++) {	//世代の全ての個体の数値の総和を求める
			total += Phenotype[i];
		}
		for (i = 0; i < IND_NUM; i++) {
			rulet_allow = (double)(Random(0, 10000)) / 10000.0 * total;	//矢の位置を決める　0~1の間で小数点４位迄の精度でランダムに小数をとり、総和をかける
			rulet_value = 0.0;
			for (j = 0; j < IND_NUM; j++) {			//矢が何番目の個体に当たったかを走査する
				rulet_value += Phenotype[j];		//番号の若い個体から順に数値を足す
				if (rulet_value > rulet_allow) {	//足しこんだ値が矢の位置を越えたら、その個体番号を記録し、次の矢を射る操作に戻る
					rulet_num = j;
					break;
				}
			}

			for (j = 0; j < GENE_LEN; j++) {	//仮の個体群配列に選択された個体を入れる
				temp_gene[i][j] = gene[rulet_num][j];
			}
		}

		for (i = 0; i < IND_NUM; i++) {			//個体数分ルーレット選択が終了したら、仮の個体を個体群に入れる
			for (j = 0; j < GENE_LEN; j++) {
				gene[i][j] = temp_gene[i][j];
			}
		}

		//交叉
		for (i = 0; i < IND_NUM; i += 2) {			//２つずつのペアで交叉する
			d_ret = (double)Random(0, 100) / 100.0;	//ペアが交叉を行うかどうかを決める　小数点以下２桁までの０～１の間の小数をランダムにとる
			if (d_ret < CROS_RATE) {				//その値が交叉率より低ければ交叉を行う
				//一点交叉
				cut = Random(0, GENE_LEN);			//一点交叉の切断点の決定（0～遺伝子長の間）
				for (j = cut; j < GENE_LEN; j++) {	//posより後の遺伝子を交叉するためforはposから開始
					temp = gene[i][j];				//遺伝子を一時的に待避
					gene[i][j] = gene[i + 1][j];	//番号の若い遺伝子を次の遺伝子に挿入
					gene[i + 1][j] = temp;			//次の遺伝子の遺伝子座に、待避させておいた遺伝子を番号の若い遺伝子に挿入
				}
			}
		}

		//突然変異
		for (i = 0; i < IND_NUM; i++) {					//各個体について全ての遺伝子座を見る
			for (j = 0; j < GENE_LEN; j++) {
				d_ret = (double)Random(0, 100) / 100.0;	//遺伝子座が突然変異を行うかどうかを決める　小数点以下２桁までの０～１の間の小数をランダムにとる
				if (d_ret < MUTA_RATE) {				//ランダムな値が突然変異率より低ければ突然変異を行う
					gene[i][j] = Random(0, ARY);		//決定した遺伝子座にARYまでのランダムな値を入れる
				}
			}
		}

		for (i = 0; i < GENE_LEN; i++) {		//エリートを最初の個体にする
			gene[0][i] = elite[i];
		}

		counter++;		//一世代に対して行う操作の終了。世代交代数をインクリメント

		//個体群の確認
		printf("第%d世代目\n", counter);
		for (i = 0; i < IND_NUM; i++) {			//個体数
			printf("%03d:", i);
			for (j = 0; j < GENE_LEN; j++) {	//遺伝子の長さ
				printf("%d", gene[i][j]);		//各遺伝子の出力
				gene_sanple[j] = gene[i][j];
			}

			printf(":%d\n", N_ary_to_decimal(gene_sanple));
		}
		printf("\n");
	}

	printf("結果:");					//最終結果の表示
	for (i = 0; i < GENE_LEN; i++) {
		printf("%d", gene[0][i]);
		gene_sanple[i] = gene[0][i];
	}
	printf(":%d\n", N_ary_to_decimal(gene_sanple));
}