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
static size_t
UKFSIF_CheckStruct(
	ukfmo_matrix_s 	*pMatrix_s_a,
	size_t 			matrixArrSize);
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
	for (j = 1u; j < (vectLen * 2u + 1u); j++)
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

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      04-сен-2019
 *
 * @brief    Функция сбрасывает массив структур указателей на матрицы
 *
 * @param[in] 	*pInit_s: 	Указатель на структуру, в которой содержится
 * 							массив указателей на матрицы
 *
 * @return  None
 */
void
UKFIS_StructInit(
	ukfsif_all_data_init_s *pInit_s)
{
	size_t i;
	for (i = 0u; i < UKFSIF_INIT_ARR_CELL_NUMB; i++)
	{
		/* Сброс параметров структуры в значения по умолчанию */
//		pInit_s->pMatrix_s_a[i]->numCols 	= 0u;
//		pInit_s->pMatrix_s_a[i]->numRows 	= 0u;
//		pInit_s->pMatrix_s_a[i]->pData 		= NULL;

		pInit_s->pMatrix_s_a[i] = NULL;
	}
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      04-сен-2019
 *
 * @brief    Функция выполняет проверку массива структур на предмет инициализации
 *
 * @para[in]  	*pMatrix_s:   	Указатель на массив структур указателей на
 * 								матрицы
 * @param[in]   matrixArrSize:	Количество ячеек в массиве структур
 * 								указателей на матрицы
 *
 * @return    Если все структуры инициализированы, то функция возвращает
 *            SIZE_MAX, иначе, возвращает индекс неинициализированной
 *            структуры в массиве
 */
static size_t
UKFSIF_CheckStruct(
	ukfmo_matrix_s 	*pMatrix_s_a,
	size_t 			matrixArrSize)
{
	size_t i;
	for (i = 0u; i < matrixArrSize; i++)
	{
		/* @FIXME добавить инкремент указателя на массив матрицы */
		if (__UKFMO_IsMatrixStructValid(pMatrix_s_a[i], SIZE_MAX, SIZE_MAX) != 1u)
		{
			return (i);
		}
	}
	return (SIZE_MAX);
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      04-сен-2019
 *
 * @brief    Функция выполняет копирование адресов структур матриц из
 *           структуры инициализации в структуру рабочих данных и проверяет
 *           на корректность
 *
 * @param[out]	*pData_s: 	Указатель на структуру рабочих данных
 * @param[in]	*pInit_s:  	Указатель на структуру инициализации
 * @param[in]  	stateLen:   Длина вектора пространства состояний
 *
 * @return  None
 */
void
UKFSIF_Init_SetMatrixPointers(
	ukfsif_all_data_s 		*pData_s,
	ukfsif_all_data_init_s 	*pInit_s,
	uint16_t 				stateLen)
{
	size_t notInitMatrixIndexNumb = SIZE_MAX;

	/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
	/* Инициализация матриц для Step 1 Calculate the sigma-points */
	pData_s->calcTheSigmaPoints_s.pMatrix_a[UKFSIF_CALC_THE_SIGMA_POINTS_x_predict] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_x_posteriori]);

	pData_s->calcTheSigmaPoints_s.pMatrix_a[UKFSIF_CALC_THE_SIGMA_POINTS_sqrtP] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_P_sqrt]);

	pData_s->calcTheSigmaPoints_s.pMatrix_a[UKFSIF_CALC_THE_SIGMA_POINTS_chi_predict] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_chi_predict]);

	pData_s->calcTheSigmaPoints_s.stateLen = stateLen;

	/* Проверка, а все ли матрицы инициализированы */
	notInitMatrixIndexNumb =
		UKFSIF_CheckStruct(
			pData_s->calcTheSigmaPoints_s.pMatrix_a[0u],
			UKFSIF_CALC_THE_SIGMA_POINTS_ARR_CELL_NUMB);
	if (notInitMatrixIndexNumb != SIZE_MAX)
	{
		/* Если попали сюда, значит одна или несколько матриц не инициализированы
		 * См. на значение notInitMatrixIndexNumb - это индекс неинициализированной структуры */
		while (1);
	}
	/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

	/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
	/* Инициализация матриц для Step 2 Calculate mean of predicted state */
	pData_s->calcMeanOfPredictState_s.meanGeneric_s.pMatrix_a[UKFSIF_CALC_MEAN_GENERIC_sigma_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_chi_apriori]);

	pData_s->calcMeanOfPredictState_s.meanGeneric_s.pMatrix_a[UKFSIF_CALC_MEAN_GENERIC_muMean] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_muMean]);

	pData_s->calcMeanOfPredictState_s.meanGeneric_s.pMatrix_a[UKFSIF_CALC_MEAN_GENERIC_vect_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_x_apriori]);

	pData_s->calcMeanOfPredictState_s.meanGeneric_s.stateLen = stateLen;

	/* Проверка, а все ли матрицы инициализированы */
	notInitMatrixIndexNumb =
		UKFSIF_CheckStruct(
			pData_s->calcMeanOfPredictState_s.meanGeneric_s.pMatrix_a[0u],
			UKFSIF_CALC_MEAN_GENERIC_ARR_CELL_NUMB);
	if (notInitMatrixIndexNumb != SIZE_MAX)
	{
		/* Если попали сюда, значит одна или несколько матриц не инициализированы
		 * См. на значение notInitMatrixIndexNumb - это индекс неинициализированной структуры */
		while (1);
	}
	/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


	/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
	/* Инициализация матриц для Step 2 Calculate covariance of predicted state  */
	pData_s->calcCovarOfPredictState_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_muCovar] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_muCovar]);

	pData_s->calcCovarOfPredictState_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_R_or_Q] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_R]);

	pData_s->calcCovarOfPredictState_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_sigma_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_chi_apriori]);

	pData_s->calcCovarOfPredictState_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_vect_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_x_apriori]);

	pData_s->calcCovarOfPredictState_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_sigma_apriori_MINUS_state_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_chi_priory_MINUS_x_priory]);

	pData_s->calcCovarOfPredictState_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_sigma_apriori_MINUS_state_apriori_TRANSPOSE] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_chi_priory_MINUS_x_priory_TRANSPOSE]);

	pData_s->calcCovarOfPredictState_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_MULT_2_MATRIX] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_result_of_mult_2_matrix]);

	pData_s->calcCovarOfPredictState_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_covariance_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_P_apriory]);

	pData_s->calcCovarOfPredictState_s.covarGeneric_s.stateLen = stateLen;

	/* Проверка, а все ли матрицы инициализированы */
	notInitMatrixIndexNumb =
		UKFSIF_CheckStruct(
			pData_s->calcCovarOfPredictState_s.covarGeneric_s.pMatrix_a[0u],
			UKFSIF_CALC_COVAR_GENERIC_ARR_CELL_NUMB);
	if (notInitMatrixIndexNumb != SIZE_MAX)
	{
		/* Если попали сюда, значит одна или несколько матриц не инициализированы
		 * См. на значение notInitMatrixIndexNumb - это индекс неинициализированной структуры */
		while (1);
	}
	/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

	/* Инициализация матриц для Step3 Propagate each sigma-point through observation
	 * Матрицы и функции, используемые на этом шаге, должны быть описаны в
	 * каждой конкретной реализации фильтра */


	/* Инициализация матриц для Step3 Calculate mean of predicted output */
	pData_s->caclMeanOfPredictOut_s.meanGeneric_s.pMatrix_a[UKFSIF_CALC_MEAN_GENERIC_sigma_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_psi_apriori]);

	pData_s->caclMeanOfPredictOut_s.meanGeneric_s.pMatrix_a[UKFSIF_CALC_MEAN_GENERIC_muMean] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_muMean]);

	pData_s->caclMeanOfPredictOut_s.meanGeneric_s.pMatrix_a[UKFSIF_CALC_MEAN_GENERIC_vect_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_y_apriori]);

	pData_s->caclMeanOfPredictOut_s.meanGeneric_s.stateLen = stateLen;

	/* Проверка, а все ли матрицы инициализированы */
	notInitMatrixIndexNumb =
		UKFSIF_CheckStruct(
			pData_s->caclMeanOfPredictOut_s.meanGeneric_s.pMatrix_a[0u],
			UKFSIF_CALC_MEAN_GENERIC_ARR_CELL_NUMB);
	if (notInitMatrixIndexNumb != SIZE_MAX)
	{
		/* Если попали сюда, значит одна или несколько матриц не инициализированы
		 * См. на значение notInitMatrixIndexNumb - это индекс неинициализированной структуры */
		while (1);
	}


	/* Инициализация матриц для Step3 Calculate covariance of predicted output */
	pData_s->caclCovarOfPredictOut_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_muCovar] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_muCovar]);

	pData_s->caclCovarOfPredictOut_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_sigma_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_psi_apriori]);

	pData_s->caclCovarOfPredictOut_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_vect_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_y_apriori]);

	pData_s->caclCovarOfPredictOut_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_sigma_apriori_MINUS_state_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_chi_priory_MINUS_x_priory]);

	pData_s->caclCovarOfPredictOut_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_sigma_apriori_MINUS_state_apriori_TRANSPOSE] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_chi_priory_MINUS_x_priory_TRANSPOSE]);

	pData_s->caclCovarOfPredictOut_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_MULT_2_MATRIX] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_result_of_mult_2_matrix]);

	pData_s->caclCovarOfPredictOut_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_covariance_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_Pyy]);

	pData_s->caclCovarOfPredictOut_s.covarGeneric_s.stateLen = stateLen;

	notInitMatrixIndexNumb =
		UKFSIF_CheckStruct(
			pData_s->caclCovarOfPredictOut_s.covarGeneric_s.pMatrix_a[0u],
			UKFSIF_CALC_COVAR_GENERIC_ARR_CELL_NUMB);
	if (notInitMatrixIndexNumb != SIZE_MAX)
	{
		/* Если попали сюда, значит одна или несколько матриц не инициализированы
		 * См. на значение notInitMatrixIndexNumb - это индекс неинициализированной структуры */
		while (1);
	}


	/* Инициализация матриц для Step3 Calculate cross-covariance of state and output */
	pData_s->calcCrossCovarOfStateAndOut_s.pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_muCovar] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_muCovar]);

	pData_s->calcCrossCovarOfStateAndOut_s.pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_chi_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_chi_apriori]);

	pData_s->calcCrossCovarOfStateAndOut_s.pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_x_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_x_apriori]);

	pData_s->calcCrossCovarOfStateAndOut_s.pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_psi_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_psi_apriori]);

	pData_s->calcCrossCovarOfStateAndOut_s.pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_y_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_y_apriori]);

	pData_s->calcCrossCovarOfStateAndOut_s.pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_chi_apriori_MINUS_x_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_chi_priory_MINUS_x_priory]);

	pData_s->calcCrossCovarOfStateAndOut_s.pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_psi_apriori_MINUS_y_apriori_TRANSPOSE] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_STEP2_psi_priory_MINUS_y_priory]);

	pData_s->calcCrossCovarOfStateAndOut_s.pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_mult_2_matrix] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_result_of_mult_2_matrix]);

	pData_s->calcCrossCovarOfStateAndOut_s.pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_Pxy] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_Pxy]);

	pData_s->calcCrossCovarOfStateAndOut_s.stateLen = stateLen;

	notInitMatrixIndexNumb =
		UKFSIF_CheckStruct(
			pData_s->calcCrossCovarOfStateAndOut_s.pMatrix_a[0u],
			UKFSIF_CALC_CROSSCOVAR_GENERIC_ARR_CELL_NUMB);
	if (notInitMatrixIndexNumb != SIZE_MAX)
	{
		/* Если попали сюда, значит одна или несколько матриц не инициализированы
		 * См. на значение notInitMatrixIndexNumb - это индекс неинициализированной структуры */
		while (1);
	}


	/* Инициализация матриц для Step4 Calculate Kalman Gain */
	pData_s->calcKalmanGain_s.pMatrix_a[UKFSIF_CALC_KALMAN_GAIN_Pxy] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_Pxy]);

	pData_s->calcKalmanGain_s.pMatrix_a[UKFSIF_CALC_KALMAN_GAIN_Pyy] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_Pyy]);

	pData_s->calcKalmanGain_s.pMatrix_a[UKFSIF_CALC_KALMAN_GAIN_Pyy_INV] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_Pyy_INV]);

	pData_s->calcKalmanGain_s.pMatrix_a[UKFSIF_CALC_KALMAN_GAIN_K] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_K]);

	pData_s->calcKalmanGain_s.stateLen = stateLen;

	notInitMatrixIndexNumb =
		UKFSIF_CheckStruct(
			pData_s->calcKalmanGain_s.pMatrix_a[0u],
			UKFSIF_CALC_KALMAN_GAIN_ARR_CELL_NUMB);
	if (notInitMatrixIndexNumb != SIZE_MAX)
	{
		/* Если попали сюда, значит одна или несколько матриц не инициализированы
		 * См. на значение notInitMatrixIndexNumb - это индекс неинициализированной структуры */
		while (1);
	}

	/* Инициализация матриц для Step4 Update state estimate */
	pData_s->updateState_s.pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_x_posteriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_x_posteriori]);

	pData_s->updateState_s.pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_x_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_x_apriori]);

	pData_s->updateState_s.pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_K] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_K]);

	pData_s->updateState_s.pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_y_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_y_apriori]);

	pData_s->updateState_s.pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_meas] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_meas]);

	pData_s->updateState_s.pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_innovation] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_innovation]);

	pData_s->updateState_s.stateLen = stateLen;

	notInitMatrixIndexNumb =
		UKFSIF_CheckStruct(
			pData_s->updateState_s.pMatrix_a[0u],
			UKFSIF_UPDATE_STATE_ESTIMATE_ARR_CELL_NUMB);
	if (notInitMatrixIndexNumb != SIZE_MAX)
	{
		/* Если попали сюда, значит одна или несколько матриц не инициализированы
		 * См. на значение notInitMatrixIndexNumb - это индекс неинициализированной структуры */
		while (1);
	}


	/* Инициализация матриц для Step4 Update error covariance */
	pData_s->updateErrCov_s.pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_P_posteriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_P]);

	pData_s->updateErrCov_s.pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_P_apriori] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_P_apriory]);

	pData_s->updateErrCov_s.pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_K] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_K]);

	pData_s->updateErrCov_s.pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_K_TRANSPOSE] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_K_TRANSPOSE]);

	pData_s->updateErrCov_s.pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_Pyy] =
		__UKFSIF_CheckMatrixStructValidation(pInit_s->pMatrix_s_a[UKFSIF_INIT_Pyy]);

	pData_s->updateErrCov_s.stateLen = stateLen;

	notInitMatrixIndexNumb =
		UKFSIF_CheckStruct(
			pData_s->updateErrCov_s.pMatrix_a[0u],
			UKFSIF_UPDATE_ERR_COVAR_ARR_CELL_NUMB);
	if (notInitMatrixIndexNumb != SIZE_MAX)
	{
		/* Если попали сюда, значит одна или несколько матриц не инициализированы
		 * См. на значение notInitMatrixIndexNumb - это индекс неинициализированной структуры */
		while (1);
	}

	/* Проверка адресов матриц */
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

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      04-сен-2019
 *
 * @brief    Функция выполняет генерацию матрицы Сигма-точек (chi_predict) на основе
 *           вектора пространства состояний (x_predict) и квадратного корня
 *           из P (cholLow(P))
 *
 * @param[in,out] 	*pData_s: 	Указатель на структуру, содержащую данные для
 * 								расчета Сигма-точек
 * @param[in]   	sqrtLAndLambda: 	Корень квадратный из (L+Lambda)
 *
 * @return  None
 */
void
UKFSIF_Step1_CalculateTheSigmaPoints(
	ukfsif_calc_the_sigma_points_s 	*pData_s,
	__UKFSIF_FPT__ 					sqrtLAndLambda)
{
	/* Заполнить первый столбец матрицы Сигма-точек */
	size_t j, i;
	for (j = 0u; (j < pData_s->stateLen); j++)
	{
		pData_s->pMatrix_a[UKFSIF_CALC_THE_SIGMA_POINTS_chi_predict]->pData[j * (pData_s->stateLen * 2u + 1u)] =
			pData_s->pMatrix_a[UKFSIF_CALC_THE_SIGMA_POINTS_x_predict]->pData[j];
	}

	UKFMO_MatrixMultScale(
		pData_s->pMatrix_a[UKFSIF_CALC_THE_SIGMA_POINTS_sqrtP],
		sqrtLAndLambda,
		pData_s->pMatrix_a[UKFSIF_CALC_THE_SIGMA_POINTS_sqrtP]);

	size_t jIdx = 0u, iIdx = 1u;
	for (j = 0u; j < pData_s->stateLen; j++)
	{
		for (i = 0u; i < pData_s->stateLen; i++)
		{
			pData_s->pMatrix_a[UKFSIF_CALC_THE_SIGMA_POINTS_chi_predict]->pData[jIdx * iIdx] =
				pData_s->pMatrix_a[UKFSIF_CALC_THE_SIGMA_POINTS_x_predict]->pData[j * i]
				+ pData_s->pMatrix_a[UKFSIF_CALC_THE_SIGMA_POINTS_sqrtP]->pData[j * i];

			pData_s->pMatrix_a[UKFSIF_CALC_THE_SIGMA_POINTS_chi_predict]->pData[jIdx * (iIdx + pData_s->stateLen)] =
				pData_s->pMatrix_a[UKFSIF_CALC_THE_SIGMA_POINTS_x_predict]->pData[j * i]
				+ pData_s->pMatrix_a[UKFSIF_CALC_THE_SIGMA_POINTS_sqrtP]->pData[j * i];

			iIdx++;
			if (iIdx >= ( pData_s->stateLen + 1u))
			{
				iIdx = 1u;
				jIdx++;
			}
		}
	}
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      04-сен-2019
 *
 * @brief    Функция выполняет расчет вектора "среднего от предсказанного состояния"
 *
 * @param[in,out] *pData_s: Указатель на структуру, содержащую данные для расчета
 *
 * @return None
 */
void
UKFSIF_Step2_CalculateMeanOfPredictedState(
	ukfsif_cacl_mean_of_predict_state_s *pData_s)
{
	UKFSIF_CaclMeanGeneric(
		&pData_s->meanGeneric_s);
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      04-сен-2019
 *
 * @brief    Функция выполняет расчет вектора "Ковариации предсказанного состояния"
 *
 * @param[in,out] *pData_s: Указатель на структуру, содержащую данные для расчета
 *
 * @return None
 */
void
UKFSIF_Step2_CalculateCovarianceOfPredictedState(
	ukfsif_calc_covar_of_predict_state_s *pData_s)
{
	/* Копирование матрицы шумов Q в матрицу "P_k|k-1" до вызова функции
	 * расчета ковариации */
	memcpy(
		(void*) pData_s->covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_covariance_apriori],
		(void*) pData_s->covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_R_or_Q],
		(size_t) pData_s->covarGeneric_s.stateLen
		* (pData_s->covarGeneric_s.stateLen * 2u + 1u));

	/* Вызов функции расчета ковариации */
	UKFSIF_CalcCovarGeneric(
		&pData_s->covarGeneric_s);
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      04-сен-2019
 *
 * @brief    Функция выполняет расчет вектора "среднего от предсказанного измерения"
 *
 * @param[in,out] *pData_s: Указатель на структуру, содержащую данные для расчета
 *
 * @return None
 */
void
UKFSIF_Step3_CalculateMeanOfPredictedOutput(
	ukfsif_calc_mean_of_predict_output_s *pData_s)
{
	UKFSIF_CaclMeanGeneric(
		&pData_s->meanGeneric_s);
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      04-сен-2019
 *
 * @brief    Функция выполняет расчет вектора "Кросс-ковариации от вектора состояния и предсказанного измерения"
 *
 * @param[in,out] *pData_s: Указатель на структуру, содержащую данные для расчета
 *
 * @return None
 */
void
UKFSIF_Step3_CalculateCrossCovarOfStateAndOut(
	ukfsif_calc_cross_covar_of_state_and_output_s *pData_s)
{
	/* Сброс матрицы кросс-ковариации в нуль */
	UKFMO_MatrixZeros(
		pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_Pxy]);

	size_t i, j;
	for (i = 0u; i < ((pData_s->stateLen * 2u) + 1u); j++)
	{
		for (j = 0u; j < pData_s->stateLen; j++)
		{
			/* Найти вектор-столбец разницы между "chi_k|k-1" и "x_k|k-1" */
			pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_chi_apriori_MINUS_x_apriori]->pData[j] =
				pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_chi_apriori]->pData[j * i]
				- pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_x_apriori]->pData[j];

			/* Найти транспонированный вектор-столбец разницы между "psi_k|k-1" и "y_k|k-1" */
			pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_psi_apriori_MINUS_y_apriori_TRANSPOSE]->pData[j] =
				pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_psi_apriori]->pData[j * i]
				- pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_y_apriori]->pData[j];
		}

		/* Умножение вектор-столбца на вектор-строку */
		UKFMO_MatrixMultiplication(
			pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_chi_apriori_MINUS_x_apriori],
			pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_psi_apriori_MINUS_y_apriori_TRANSPOSE],
			pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_mult_2_matrix]);

		/* Умножение матрицы на скалярный коэффициент из вектора весовых
		 * коэффициентов */
		UKFMO_MatrixMultScale(
			pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_mult_2_matrix],
			pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_muCovar]->pData[i],
			pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_mult_2_matrix]);

		/* Сложить с результатом, полученным на предыдущей итерации цикла */
		UKMO_MatrixAdition(
			pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_Pxy],
			pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_mult_2_matrix],
			pData_s->pMatrix_a[UKFSIF_CALC_CROSSCOVAR_GENERIC_Pxy]);
	}
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      04-сен-2019
 *
 * @brief    Функция выполняет расчет вектора "Ковариация от предсказанного измерения"
 *
 * @param[in,out] *pData_s: Указатель на структуру, содержащую данные для расчета
 *
 * @return None
 */
void
UKFSIF_Step3_CalculateCovarianceOfPredictedOutput(
	ukfsif_calc_covar_of_predict_output_s *pData_s)
{
	UKFMO_CopyMatrix(
		pData_s->covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_covariance_apriori],
		pData_s->covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_R_or_Q]);

	UKFSIF_CalcCovarGeneric(
		&pData_s->covarGeneric_s);
}

void
UKFSIF_Step3_CalculateCovarianceOfPredictedOutputAndCrossCovariance(
	ukfsif_all_data_s *pData_s)
{
	/* Calculate covariance of predicted output -->>>>>>>>>>>>>>>>>>>>>>>>>>> */
	/* Скопировать матрицу шумов в матрицу Ковариации */
	UKFMO_CopyMatrix(
		pData_s->caclCovarOfPredictOut_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_covariance_apriori],
		pData_s->caclCovarOfPredictOut_s.covarGeneric_s.pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_R_or_Q]);

	/* Рассчитать ковариацию */
	UKFSIF_CalcCovarGeneric(
		&pData_s->caclCovarOfPredictOut_s.covarGeneric_s);
	/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


	/* Calculate cross-covariance of state and output -->>>>>>>>>>>>>>>>>>>>> */
	UKFSIF_Step3_CalculateCrossCovarOfStateAndOut(
		&pData_s->calcCrossCovarOfStateAndOut_s);

	/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      04-сен-2019
 *
 * @brief    Универсальная функция для расчета Ковариации
 *
 * @param[in,out] *pData_s: Указатель на структуру, содержащую данные для расчета
 *
 * @return None
 */
void
UKFSIF_CalcCovarGeneric(
	ukfsif_calc_covar_generic_s *pData_s)
{
	size_t i, j;
	for (i = 0u; i < ((pData_s->stateLen * 2u) + 1u); j++)
	{
		/* Найти разницу между вектор-столбцом пространства состояний и
		 * вектор-столбцом матрицы Сигма-точек */
		for (j = 0u; j < pData_s->stateLen; j++)
		{
			pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_sigma_apriori_MINUS_state_apriori]->pData[j] =
				pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_sigma_apriori]->pData[j * i]
				- pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_vect_apriori]->pData[j * pData_s->stateLen];

			/* Копирование каждого столбца во временный массив
			 * "UKFSIF_STEP2_CHI_MINUS_STATE_TEMP" размерностью (Lx2L+1)
			 * для дальнейшего использования на "Step3 Calculate cross-covariance
			 * of state and output" */
			// pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_sigma_apriori_MINUS_state_apriori_TEMP]->pData[i * j] =
			// 	pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_sigma_apriori_MINUS_state_apriori]->pData[j];

		}

		/* Транспонирование */
		UKFMO_MatrixTranspose(
			pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_sigma_apriori_MINUS_state_apriori],
			pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_sigma_apriori_MINUS_state_apriori_TRANSPOSE]);

		/* Умножение вектор-столбца на его транспонированную версию */
		UKFMO_MatrixMultiplication(
			pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_sigma_apriori_MINUS_state_apriori],
			pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_sigma_apriori_MINUS_state_apriori_TRANSPOSE],
			pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_MULT_2_MATRIX]);

		/* Умножение матрицы LxL на скаляр весового коэффициента */
		UKFMO_MatrixMultScale(
			pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_MULT_2_MATRIX],
			pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_muCovar]->pData[i],
			pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_MULT_2_MATRIX]);

		/* Сложить с предыдущим результатом */
		UKMO_MatrixAdition(
			pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_covariance_apriori],
			pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_MULT_2_MATRIX],
			pData_s->pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_covariance_apriori]);
	}
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      04-сен-2019
 *
 * @brief    Универсальная функция для расчета "Среднего вектора"
 *
 * @param[in,out] *pData_s: Указатель на структуру, содержащую данные для расчета
 *
 * @return  None
 */
void
UKFSIF_CaclMeanGeneric(
	ukfsif_calc_mean_generic_s *pData_s)
{
	UKFMO_MatrixMultiplication(
		pData_s->pMatrix_a[UKFSIF_CALC_MEAN_GENERIC_sigma_apriori],
		pData_s->pMatrix_a[UKFSIF_CALC_MEAN_GENERIC_muMean],
		pData_s->pMatrix_a[UKFSIF_CALC_MEAN_GENERIC_vect_apriori]);
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      04-сен-2019
 *
 * @brief    Функция выполняет расчет матрицы усиления
 *
 * @param[in,out] *pData_s: Указатель на структуру, содержащую данные для расчета
 *
 * @return  None
 */
void
UKFSIF_Step4_CalcKalmanGain(
	ukfsif_calc_kalman_gain_s *pData_s)
{
	/* Найти обратную матрицу от P_yy */
	UKFMO_MatrixIdentity(
		pData_s->pMatrix_a[UKFSIF_CALC_KALMAN_GAIN_Pyy_INV]);
	UKFMO_MatrixInverse(
		pData_s->pMatrix_a[UKFSIF_CALC_KALMAN_GAIN_Pyy],
		pData_s->pMatrix_a[UKFSIF_CALC_KALMAN_GAIN_Pyy_INV]);

	/* Найти матрицу усиления Калмана */
	UKFMO_MatrixMultiplication(
		pData_s->pMatrix_a[UKFSIF_CALC_KALMAN_GAIN_Pxy],
		pData_s->pMatrix_a[UKFSIF_CALC_KALMAN_GAIN_Pyy_INV],
		pData_s->pMatrix_a[UKFSIF_CALC_KALMAN_GAIN_K]);

}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      04-сен-2019
 *
 * @brief    Функция выполняет коррекцию вектора пространства состояний с
 *           помощью внешнего измерителя
 *
 * @param[in,out] *pData_s: Указатель на структуру, содержащую данные для
 * 							расчета
 *
 * @return  None
 */
void
UKFSIF_Step4_UpdateStateEstimate(
	ukfsif_update_state_s *pData_s)
{
	/* Найти Инновацию, т.е. разницу между вектором измерений и его
	 * предсказанным значением */
	UKMO_MatrixSubstraction(
		pData_s->pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_meas],
		pData_s->pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_x_apriori],
		pData_s->pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_innovation]);

	/* Умножить матрицу усиления на Инновацию */
	UKFMO_MatrixMultiplication(
		pData_s->pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_K],
		pData_s->pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_innovation],
		pData_s->pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_innovation]);

	/* Сложить вектор пространства состояний с умноженной на матрицу
	 * усиления Инновацией */
	UKMO_MatrixAdition(
		pData_s->pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_x_apriori],
		pData_s->pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_innovation],
		pData_s->pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_x_posteriori]);
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      04-сен-2019
 *
 * @brief    Функция обновляет значение матрицы ковариации (P)
 *
 * @param[in,out] *pData_s: Указатель на структуру, содержащую данные для
 * 							расчета
 *
 * @return  None
 */
void
UKFSIF_Step4_UpdateErrorCovariance(
	ukfsif_update_err_covar_s *pData_s)
{
	/* Найти транспонированную матрицу от матрицы коэффициентов усиления */
	UKFMO_MatrixTranspose(
		pData_s->pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_K],
		pData_s->pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_K_TRANSPOSE]);

	/* Умножить матрицу усиления на Pyy */
	UKFMO_MatrixMultiplication(
		pData_s->pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_K],
		pData_s->pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_Pyy],
		pData_s->pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_P_posteriori]);

	/* Результат умножения K и Pyy  умножить на транспонированную матрицу
	 * коэффициентов усиления */
	UKFMO_MatrixMultiplication(
		pData_s->pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_P_posteriori],
		pData_s->pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_K_TRANSPOSE],
		pData_s->pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_P_posteriori]);

	/* из "P_k|k-1" вычесть результат, полученный выше */
	UKMO_MatrixSubstraction(
		pData_s->pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_P_apriori],
		pData_s->pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_P_posteriori],
		pData_s->pMatrix_a[UKFSIF_UPDATE_ERR_COVAR_P_posteriori]);

}
/*#### |End  | <-- Секция - "Описание глобальных функций" ####################*/


/*#### |Begin| --> Секция - "Создание задач" #################################*/
/*#### |End  | <-- Секция - "Создание задач" #################################*/


/*#### |Begin| --> Секция - "Описание локальных функций" #####################*/
/*#### |End  | <-- Секция - "Описание локальных функций" #####################*/


/*############################################################################*/
/*############################ END OF FILE  ##################################*/
/*############################################################################*/
