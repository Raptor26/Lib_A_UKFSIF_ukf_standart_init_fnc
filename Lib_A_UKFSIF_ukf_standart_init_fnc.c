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
/*#### |End  | <-- Секция - "Описание глобальных функций" ####################*/


/*#### |Begin| --> Секция - "Создание задач" #################################*/
/*#### |End  | <-- Секция - "Создание задач" #################################*/


/*#### |Begin| --> Секция - "Описание локальных функций" #####################*/
/*#### |End  | <-- Секция - "Описание локальных функций" #####################*/


/*############################################################################*/
/*############################ END OF FILE  ##################################*/
/*############################################################################*/
