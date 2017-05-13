#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//��̐�����\������̐� (�����̏�����A�����ɂ��邱��)
#define IND_NUM 20
//��`�q�̒����@��������ƌv�Z���ɃI�[�o�[�t���[����̂Œ���
#define GENE_LEN 16
//N�i�@�i�����ɂ͎w�肵�����i�@����1���������ɂ��邱�Ɓ@�傫������ƌv�Z���ɃI�[�o�[�t���[����̂Œ��Ӂj
#define ARY 1
//������̉�
#define SEDAI_MAX 50
//������
#define CROS_RATE 0.95
//�ˑR�ψٗ�
#define MUTA_RATE 0.05

int Random(int min, int max)			//���������֐� �����m���󂯎��O����(�m�|�P)�܂ł̐����������_���ɕԂ�
{
	return min + (int)(rand()*(max - min + 1.0) / (1.0 + RAND_MAX));
}

int N_ary_to_decimal(int gene_sanmple[GENE_LEN])		//N�i�@����10�i�@�ɕϊ�����֐�
{
	int decimal = 0;	//10�i�@
	int digit;			//����
	int base = 1;		//��

	for (digit = GENE_LEN; digit >0; digit--) {

		decimal = decimal + gene_sanmple[digit - 1] * base;
		base = base * (ARY + 1);
	}
	return decimal;
}

void main() {
	int i, j;						//�v�Z�p�ϐ�
	int counter;					//�����㐔

	int gene[IND_NUM][GENE_LEN];	//�̌Q�̔z��
	int gene_sanple[GENE_LEN];		//1�̂̈�`�q��i�v�Z�p�j
	int Phenotype[IND_NUM];			//�̂̐��l

	int elite_Ind;					//�G���[�g�̂̔ԍ�
	int elite_num;					//�G���[�g�̂̐��l
	int elite[GENE_LEN];			//�G���[�g�̂̕ۑ��p�z��

	double total, rulet_value, rulet_allow;	//���[���b�g�I��p�ϐ�
	int rulet_num;							//���[���b�g�őI�����ꂽ�̂̔ԍ�
	int temp_gene[IND_NUM][GENE_LEN];		//���[���b�g�I���ŗp����̌Q�̑ޔ�p�z��

	int cut;						//�����ɂ�����ؒf�_��\���ϐ�
	int temp;						//�����ɂ�������֗p�ϐ�
	double d_ret;					//�����y�ѓˑR�ψقŗp����ϐ�

	time_t t;						//���������p�̕ϐ�

	srand((unsigned)time(&t));	//time�֐��𗘗p���ė����̐���

	//�̌Q�̐���
	for (i = 0; i < IND_NUM; i++) {			//�̐�
		for (j = 0; j < GENE_LEN; j++) {	//��`�q�̒���
			gene[i][j] = Random(0, ARY);	//�e��`�q�Ƀ����_���Ɉ�`�q�����
		}

	}
	
	//�������ꂽ�̌Q�̊m�F
	printf("�������ꂽ�̌Q\n");
	for (i = 0; i < IND_NUM; i++) {			//�̐�
		printf("%03d:", i);
		for (j = 0; j < GENE_LEN; j++) {	//��`�q�̒���
			printf("%d", gene[i][j]);		//�e��`�q�̏o��
			gene_sanple[j] = gene[i][j];
		}
		printf(":%d\n", N_ary_to_decimal(gene_sanple));
	}
	printf("\n");
	
	counter = 0;	//�����㐔�����Z�b�g

	while (counter < SEDAI_MAX) {		//�����㐔�ɒB����܂Ń��[�v

		//�̂̕]���@�����10�i���ɂ���Ƃ��傫�����ǂ��Ƃ��������Ƃ���B
		for (i = 0; i < IND_NUM; i++) {			//�̐������[�v
			for (j = 0; j < GENE_LEN; j++) {	//��`�q�̒��������[�v
				gene_sanple[j] = gene[i][j];	//�̌Q��1�̈�`�q��z��Ɋi�[
			}	
			
			Phenotype[i] = N_ary_to_decimal(gene_sanple);	//�e�̂̐��l���v�Z
		}

		//�G���[�g�̂̑���
		elite_Ind = 0;		//�ŏ����G���[�g�Ɖ���
		elite_num = Phenotype[0];			//��ԍŏ��̌̂̐��т���
		for (i = 1; i < IND_NUM; i++) {		//�S�Ă̌̂̐��l�ƃG���[�g�̐��l���ׂ�
			if (elite_num < Phenotype[i]) {	//���̃G���[�g�������т̂悢�̂�����������G���[�g�̂̔ԍ��Ɛ��l���X�V����
				elite_Ind = i;
				elite_num = Phenotype[i];
			}
		}
		for (i = 0; i < GENE_LEN; i++) {	//elite_Ind�̈�`�q���G���[�g�z��ɓ����
			elite[i] = gene[elite_Ind][i];
		}


		//���[���b�g�I��
		total = 0.0;					//���a�̃��Z�b�g
		for (i = 0; i < IND_NUM; i++) {	//����̑S�Ă̌̂̐��l�̑��a�����߂�
			total += Phenotype[i];
		}
		for (i = 0; i < IND_NUM; i++) {
			rulet_allow = (double)(Random(0, 10000)) / 10000.0 * total;	//��̈ʒu�����߂�@0~1�̊Ԃŏ����_�S�ʖ��̐��x�Ń����_���ɏ������Ƃ�A���a��������
			rulet_value = 0.0;
			for (j = 0; j < IND_NUM; j++) {			//����Ԗڂ̌̂ɓ����������𑖍�����
				rulet_value += Phenotype[j];		//�ԍ��̎Ⴂ�̂��珇�ɐ��l�𑫂�
				if (rulet_value > rulet_allow) {	//�������񂾒l����̈ʒu���z������A���̌̔ԍ����L�^���A���̖���˂鑀��ɖ߂�
					rulet_num = j;
					break;
				}
			}

			for (j = 0; j < GENE_LEN; j++) {	//���̌̌Q�z��ɑI�����ꂽ�̂�����
				temp_gene[i][j] = gene[rulet_num][j];
			}
		}

		for (i = 0; i < IND_NUM; i++) {			//�̐������[���b�g�I�����I��������A���̌̂��̌Q�ɓ����
			for (j = 0; j < GENE_LEN; j++) {
				gene[i][j] = temp_gene[i][j];
			}
		}

		//����
		for (i = 0; i < IND_NUM; i += 2) {			//�Q���̃y�A�Ō�������
			d_ret = (double)Random(0, 100) / 100.0;	//�y�A���������s�����ǂ��������߂�@�����_�ȉ��Q���܂ł̂O�`�P�̊Ԃ̏����������_���ɂƂ�
			if (d_ret < CROS_RATE) {				//���̒l�����������Ⴏ��Ό������s��
				//��_����
				cut = Random(0, GENE_LEN);			//��_�����̐ؒf�_�̌���i0�`��`�q���̊ԁj
				for (j = cut; j < GENE_LEN; j++) {	//pos����̈�`�q���������邽��for��pos����J�n
					temp = gene[i][j];				//��`�q���ꎞ�I�ɑҔ�
					gene[i][j] = gene[i + 1][j];	//�ԍ��̎Ⴂ��`�q�����̈�`�q�ɑ}��
					gene[i + 1][j] = temp;			//���̈�`�q�̈�`�q���ɁA�Ҕ������Ă�������`�q��ԍ��̎Ⴂ��`�q�ɑ}��
				}
			}
		}

		//�ˑR�ψ�
		for (i = 0; i < IND_NUM; i++) {					//�e�̂ɂ��đS�Ă̈�`�q��������
			for (j = 0; j < GENE_LEN; j++) {
				d_ret = (double)Random(0, 100) / 100.0;	//��`�q�����ˑR�ψق��s�����ǂ��������߂�@�����_�ȉ��Q���܂ł̂O�`�P�̊Ԃ̏����������_���ɂƂ�
				if (d_ret < MUTA_RATE) {				//�����_���Ȓl���ˑR�ψٗ����Ⴏ��ΓˑR�ψق��s��
					gene[i][j] = Random(0, ARY);		//���肵����`�q����ARY�܂ł̃����_���Ȓl������
				}
			}
		}

		for (i = 0; i < GENE_LEN; i++) {		//�G���[�g���ŏ��̌̂ɂ���
			gene[0][i] = elite[i];
		}

		counter++;		//�ꐢ��ɑ΂��čs������̏I���B�����㐔���C���N�������g

		//�̌Q�̊m�F
		printf("��%d�����\n", counter);
		for (i = 0; i < IND_NUM; i++) {			//�̐�
			printf("%03d:", i);
			for (j = 0; j < GENE_LEN; j++) {	//��`�q�̒���
				printf("%d", gene[i][j]);		//�e��`�q�̏o��
				gene_sanple[j] = gene[i][j];
			}

			printf(":%d\n", N_ary_to_decimal(gene_sanple));
		}
		printf("\n");
	}

	printf("����:");					//�ŏI���ʂ̕\��
	for (i = 0; i < GENE_LEN; i++) {
		printf("%d", gene[0][i]);
		gene_sanple[i] = gene[0][i];
	}
	printf(":%d\n", N_ary_to_decimal(gene_sanple));
}