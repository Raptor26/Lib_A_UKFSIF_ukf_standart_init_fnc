/**
 * @file   	%<%NAME%>%.%<%EXTENSION%>%
 * @author 	%<%USER%>%
 * @version
 * @date 	%<%DATE%>%, %<%TIME%>%
 * @brief
 */


/*#### |Begin| --> Секция - "Include" ########################################*/
#include "Lib_A_UKFSIF_ukf_standart_init_fnc.h"
/*#### |End  | <-- Секция - "Include" ########################################*/


/*#### |Begin| --> Секция - "Глобальные переменные" ##########################*/
/*#### |End  | <-- Секция - "Глобальные переменные" ##########################*/


/*#### |Begin| --> Секция - "Локальные переменные" ###########################*/
/*#### |End  | <-- Секция - "Локальные переменные" ###########################*/


/*#### |Begin| --> Секция - "Прототипы локальных функций" ####################*/
/*#### |End  | <-- Секция - "Прототипы локальных функций" ####################*/


/*#### |Begin| --> Секция - "Описание глобальных функций" ####################*/

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      22-авг-2019
 *
 * @brief    Функция выполняет расчет скалярного параметра "Lambda"
 *
 * @param[in]    vectLen: 	Длина вектора пространства состояний
 * @param[in]    alpha:     Коэффициент
 * @param[in]    kappa:     Коэффициент
 *
 * @return Рассчитанный скалярный параметр "Lambda"
 */
__UKFSIF_FPT__
UKFSIF_GetLambda(
	uint16_t 		vectLen,
	__UKFSIF_FPT__ 	alpha,
	__UKFSIF_FPT__ 	kappa)
{
	return (((alpha *  alpha) * (((__UKFSIF_FPT__) vectLen) + kappa)) - ((__UKFSIF_FPT__) vectLen));
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      22-авг-2019
 *
 * @brief    Функция выполняет расчет вектора весовых коэффициентов "Mean"
 *
 * @param[in] 	*pScalParams_s:	Указатель на структуру, содержащую
 * 								скалярные параметры
 * @param[out] 	*pWeightMean: 	Указатель на область памяти, в которую
 * 								будет записан вектор весовых коэффициентов
 * @param[in]    vectLen: 	Длина вектора пространства состояний
 *
 * @return  None
 */
void
UKFSIF_InitWeightVectorMean(
	ukfsif_scaling_param_s 	*pScalParams_s,
	__UKFSIF_FPT__ 			*pWeightMean,
	uint16_t 				 vectLen)
{
	/* Calculate scaling parameter */
	__UKFSIF_FPT__ lambda =
		UKFSIF_GetLambda(
			vectLen,
			pScalParams_s->alpha,
			pScalParams_s->kappa);

	/* Calculate weight vector Mean */
	size_t j;
	__UKFSIF_FPT__ eta_meanAndCovar =
		((__UKFSIF_FPT__) 1.0) / (((__UKFSIF_FPT__)2.0) * (((__UKFSIF_FPT__) vectLen) + lambda));
	for (j = 1u; j < vectLen * 2u; j++)
	{
		pWeightMean[j] = eta_meanAndCovar;
	}

	pWeightMean[0u] =
		lambda / (((__UKFSIF_FPT__) vectLen) + lambda);
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      22-авг-2019
 *
 * @brief    Функция выполняет расчет вектора весовых коэффициентов "Cov"
 *
 * @param[in] 	*pScalParams_s:	Указатель на структуру, содержащую
 * 								скалярные параметры
 * @param[out] 	*pWeightCov: 	Указатель на область памяти, в которую
 * 								будет записан вектор весовых коэффициентов
 * @param[in]    vectLen: 	Длина вектора пространства состояний
 *
 * @return  None
 */
void
UKFSIF_InitWeightVectorCov(
	ukfsif_scaling_param_s 	*pScalParams_s,
	__UKFSIF_FPT__ 			*pWeightCov,
	uint16_t 				 vectLen)
{
	/* Формулы для расчета весовых коэффициентов "Mean" и "Cov" совпадают
	 * кроме нулевых элементов вектора, поэтому для расчета вектора "Cov"
	 * используется функция, которая считает "Mean", а затем пересчитывается
	 * только нулевой элемент вектора "Cov" */

	/* Заполнить значениями как и для вектора "Mean" */
	UKFSIF_InitWeightVectorMean(
		pScalParams_s,
		pWeightCov,
		vectLen);

	/* Пересчитать только нулевой элемент вектора "Cov": */
	/* -найти Lambda */
	__UKFSIF_FPT__ lambda =
		UKFSIF_GetLambda(
			vectLen,
			pScalParams_s->alpha,
			pScalParams_s->kappa);

	/* -и посчитать нулевой элемент вектора "Cov" */
	pWeightCov[0u] =
		(lambda / (((__UKFSIF_FPT__) vectLen) + lambda)) + ((__UKFSIF_FPT__)1.0) - (pScalParams_s->alpha * pScalParams_s->alpha) + pScalParams_s->beta;
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      26-авг-2019
 *
 * @brief    Функция генерирует массив Сигма-точек методом 2L+1
 *
 * @param[in]       *pStateVect:	Указатель на двумерный массив,
 * 									содержащий матрицу пространства состояний
 * @param[out]   	*pSigmaPoints: Указатель на двумерный массив,
 * @param[in,out]   *pSqrtP: Пп sqrt p
 * @param[in]    	sqrtLenLambda: Длина лямбда
 * @param[in]    	stateVectLen: Государство vect длина
 * @param[in]    	sigmaPointsColNumb:    Сигма указывает на онемение
 */
void
UKFSIF_CalculateTheSigmaPoints_2L1(
	__UKFSIF_FPT__ *pStateVect,
	__UKFSIF_FPT__ *pSigmaPoints,
	__UKFSIF_FPT__ *pSqrtP, 			/* Указатель на двумерный массив, в котором содержится квадратный корень из матрицы ковариации */
	__UKFSIF_FPT__  sqrtLenLambda,
	uint16_t 		stateVectLen 		/* Длина вектора пространства состояний, совпадает с количеством строк матрицы Сигма-точек */
)
{
	/* Количество столбцов матрицы Сигма-точек */
	uint16_t sigmaPointsColNumb =
		(stateVectLen * 2u) + 1u;

	/* Заполнение 1-го столбца матрицы Сигма-точек */
	size_t j;
	for (j = 0u; j < stateVectLen; j++)
	{
		/* Копирование матрицы вектора пространства состояний в 1-й
		 * столбец матрицы Сигма-точек */
		/* Помни, pSigmaPoints 	- указатель на двумерный массив
		 * Помни, pStateVect 	- указатель на двумерный массив */
		pSigmaPoints[j * stateVectLen] =
			pStateVect[j * stateVectLen];
	}

	/* Умножить матрицу квадратного корня из ковариации на скаляр */
	for (j = 0u; j < (stateVectLen * sigmaPointsColNumb); j++)
	{
		pSqrtP[j] *= sqrtLenLambda;
	}

	/* Генерация остальных сигма-точек */
	size_t i;
	size_t jIdx = 0u, iIdx = 1u;
	/* @FIXME Удалить эту переменную и заменить на значение из принимаемого параметра */
	size_t iIdxOffset = stateVectLen;
	for (j = 0u; j < stateVectLen; j++)
	{
		for (i = 0u; i < stateVectLen; i++)
		{
			pSigmaPoints[jIdx * iIdx] =
				pStateVect[j * i] + pSqrtP[j * i];

			pSigmaPoints[jIdx * (iIdx + iIdxOffset)] =
				pStateVect[j * i] - pSqrtP[j * i];
			iIdx++;
			if (iIdx >= (iIdxOffset + 1u))
			{
				iIdx = 1u;
				jIdx++;
			}
		}
	}
}

void
UKFSIF_StructInit_Step2Data(
	ukfsif_step2_params_2l1_s 		*pData_s)
{
	size_t i;
	for (i = 0; i < UKFSIF_STEP2_ARR_CELL_NUMB; i++)
	{
		/* Сброс указателей на структуры в NULL */
		pData_s->pMatrix_a[i] = NULL;
	}

	/* Сброс длины пространства состояний в нуль */
	pData_s->stateLen = 0u;
}

void
UKFSIF_Init_Step2Data(
	ukfsif_step2_params_2l1_s 		*pData_s,
	ukfsif_step2_params_2l1_init_s 	*pInit_s)
{
	/* Адреса структур "ukfsif_step2_params_2l1_s" и
	 * "ukfsif_step2_params_2l1_init_s" не должны совпадать */
	if (pInit_s == pData_s)
	{
		while (1);
	}

//	/* Сброс указателей в нуль */
//	size_t i;
//	for (i = 0; i < UKFSIF_STEP2_MEM_CELL_NUMB; i++)
//	{
//		pData_s->pMemAddr_a[i] = NULL;
//	}
//	/* Сброс длины вектора пространства состояний */
//	pData_s->stateLen = 0u;
//
//	/* Копирование указателей из структуры инициализации в рабочую структуру */
//	for (i = 0; i < UKFSIF_STEP2_MEM_CELL_NUMB; i++)
//	{
//		pData_s->pMemAddr_a[i] = pInit_s->pMemAddr_a[i];
//	}
//
//	/* Копирование длины вектора пространства состояний */
//	pData_s->stateLen = pInit_s->stateLen;
//
//	/* Проверка области памяти на указатель NULL */
//	for (i = 0; i < UKFSIF_STEP2_MEM_CELL_NUMB; i++)
//	{
//		/* Если адрес равен NULL */
//		if (pData_s->pMemAddr_a[i] == NULL)
//		{
//			/* Зависнуть */
//			while (1);
//		}
//	}
//	/* Если длина вектора пространства состояний меньше 1 */
//	if (pData_s->stateLen < 1u)
//	{
//		/* Зависнуть */
//		while (1);
//	}

	size_t i;

	/* Проверка массива указателей на структуры */
	for (i = 0;
		 i < ((size_t) UKFSIF_STEP2_ARR_CELL_NUMB);
		 i++)
	{
		/* Если один из указателей на матрицы не инициализирован */
		if (pData_s->pMatrix_a[i]->pData 	== NULL ||
			pData_s->pMatrix_a[i]->numCols 	== 0u 	||
			pData_s->pMatrix_a[i]->numRows 	== 0u	||
			pData_s->pMatrix_a[i] 			== NULL)
		{
			/* Зависнуть */
			while (1);
		}
	}
}

void
UKFSIF_Step2_CalculateMeanOfPredictedState(
	ukfsif_step2_params_2l1_s *pData_s)
{
	UKFMO_MatrixMultiplication(
		pData_s->pMatrix_a[UKFSIF_STEP2_chi_priory],
		pData_s->pMatrix_a[UKFSIF_STEP2_MUMEAN],
		pData_s->pMatrix_a[UKFSIF_STEP2_x_priory]);
}

void
UKFSIF_Step2_CalculateCovarianceOfPredictedState(
	ukfsif_step2_params_2l1_s *pData_s)
{
	/* Копирование матрицы шумов в матрицу Ковариации "P_k|k-1" */
	memcpy(
		(void*) pData_s->pMatrix_a[UKFSIF_STEP2_P_apriory]->pData,
		(void*) pData_s->pMatrix_a[UKFSIF_STEP2_Q]->pData,
		pData_s->stateLen * ((pData_s->stateLen * 2) + 1u));

	size_t i, j;
	for (i = 0u; i < ((pData_s->stateLen * 2u) + 1u); j++)
	{
		/* Найти разницу между вектор-столбцом пространства состояний и
		 * вектор-столбцом матрицы Сигма-точек */
		for (j = 0u; j < pData_s->stateLen; j++)
		{
			pData_s->pMatrix_a[UKFSIF_STEP2_chi_priory_MINUS_x_priory]->pData[j] =
				pData_s->pMatrix_a[UKFSIF_STEP2_chi_priory]->pData[j * i]
				- pData_s->pMatrix_a[UKFSIF_STEP2_x_priory]->pData[j * pData_s->stateLen];

			/* Копирование каждого столбца во временный массив
			 * "UKFSIF_STEP2_CHI_MINUS_STATE_TEMP" размерностью (Lx2L+1)
			 * для дальнейшего использования на "Step3 Calculate cross-covariance
			 * of state and output" */
			pData_s->pMatrix_a[UKFSIF_STEP2_chi_priory_MINUS_x_priory_TEMP]->pData[i * j] =
				pData_s->pMatrix_a[UKFSIF_STEP2_chi_priory_MINUS_x_priory]->pData[j];
		}

		/* Транспонирование */
		UKFMO_MatrixTranspose(
			pData_s->pMatrix_a[UKFSIF_STEP2_chi_priory_MINUS_x_priory],
			pData_s->pMatrix_a[UKFSIF_STEP2_chi_priory_MINUS_x_priory_TRANPOSE]);

		/* Умножение вектор-столбца на его транспонированную версию */
		UKFMO_MatrixMultiplication(
			pData_s->pMatrix_a[UKFSIF_STEP2_chi_priory_MINUS_x_priory],
			pData_s->pMatrix_a[UKFSIF_STEP2_chi_priory_MINUS_x_priory_TRANPOSE],
			pData_s->pMatrix_a[UKFSIF_STEP2_RESULT_OF_MULT_2_MATRIX]);

		/* Умножение матрицы LxL на скаляр весового коэффициента */
		UKFMO_MatrixMultScale(
			pData_s->pMatrix_a[UKFSIF_STEP2_RESULT_OF_MULT_2_MATRIX],
			pData_s->pMatrix_a[UKFSIF_STEP2_MUCOV]->pData[i],
			pData_s->pMatrix_a[UKFSIF_STEP2_RESULT_OF_MULT_2_MATRIX]);

		/* Сложить с предыдущим результатом */
		UKMO_MatrixAdition(
			pData_s->pMatrix_a[UKFSIF_STEP2_P_apriory],
			pData_s->pMatrix_a[UKFSIF_STEP2_RESULT_OF_MULT_2_MATRIX],
			pData_s->pMatrix_a[UKFSIF_STEP2_P_apriory]);
	}

}
/*#### |End  | <-- Секция - "Описание глобальных функций" ####################*/


/*#### |Begin| --> Секция - "Создание задач" #################################*/
/*#### |End  | <-- Секция - "Создание задач" #################################*/


/*#### |Begin| --> Секция - "Описание локальных функций" #####################*/
/*#### |End  | <-- Секция - "Описание локальных функций" #####################*/


/*############################################################################*/
/*############################ END OF FILE  ##################################*/
/*############################################################################*/
