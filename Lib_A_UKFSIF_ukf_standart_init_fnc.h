/**
 * @file   	%<%NAME%>%.%<%EXTENSION%>%
 * @author 	%<%USER%>%
 * @version
 * @date 	%<%DATE%>%, %<%TIME%>%
 * @brief 	Правила присваивания имен матрицам:
 *			"k-1"	- "predict"
 *			"k|k-1"	- "priory"
 *			"k"		- "postreory"
 */


#ifndef LIB_A_UKFSIF_UKF_STANDART_INIT_FNC_H_
#define LIB_A_UKFSIF_UKF_STANDART_INIT_FNC_H_


/*#### |Begin| --> Секция - "Include" ########################################*/
/*==== |Begin| --> Секция - "C libraries" ====================================*/
#include "stdint.h"
#include "stdio.h"
/*==== |End  | <-- Секция - "C libraries" ====================================*/

/*==== |Begin| --> Секция - "RTOS libraries ==================================*/
/*==== |End  | <-- Секция - "RTOS libraries ==================================*/

/*==== |Begin| --> Секция - "MK peripheral libraries" ========================*/
/*==== |End  | <-- Секция - "MK peripheral libraries" ========================*/

/*==== |Begin| --> Секция - "Extern libraries" ===============================*/
#include "Lib_A_UKFMO_ukf_matrix_operations.h"
/*==== |End  | <-- Секция - "Extern libraries" ===============================*/
/*#### |End  | <-- Секция - "Include" ########################################*/


/*#### |Begin| --> Секция - "Определение констант" ###########################*/

/*==== |Begin| --> Секция определения типа числа с плавающей точкой ==========*/
#if !defined (__UKFSIF_FPT__)
	#error "Please, set __UKFSIF_FPT__ float or double in macros list"
#endif

#if !defined (__UKFSIF_FPT_SIZE__)
	#error "Please, set __UKFSIF_FPT_SIZE__ 4 (that mean float) or 8 (that mean double) in macros list"
#endif

#if     __UKFSIF_FPT_SIZE__ == 4

#elif   __UKFSIF_FPT_SIZE__ == 8

#else
	#error "Your compiler uses a non-standard floating point size"
#endif
/*==== |End  | <-- Секция определения типа числа с плавающей точкой ==========*/

/*==== |Begin| --> Секция - Макросы для встраиваемых функций =================*/
#if defined (__GNUC__)

	/* inline*/
	#ifndef __UKFSIF_INLINE
		#define __UKFSIF_INLINE          	inline
	#endif

	/* static inline */
	#ifndef __UKFSIF_STATIC_INLINE
		#define __UKFSIF_STATIC_INLINE   	static inline
	#endif

	/* always inline */
	#ifndef __UKFSIF_ALWAYS_INLINE
		#define __UKFSIF_ALWAYS_INLINE    	inline __attribute__((always_inline)) static
	#endif

	/* force inline */
	#ifndef __UKFSIF_FORCE_INLINE
		#define __UKFSIF_FORCE_INLINE    	inline __attribute__((always_inline))
	#endif

#else
	#define __UKFSIF_INLINE
	#define __UKFSIF_STATIC_INLINE   static
	#define __UKFSIF_ALWAYS_INLINE
#endif
/*==== |End  | <-- Секция - Макросы для встраиваемых функций =================*/

/*#### |End  | <-- Секция - "Определение констант" ###########################*/


/*#### |Begin| --> Секция - "Определение типов" ##############################*/

typedef enum
{
	UKFSIF_STEP2_MUMEAN = 0u,
	UKFSIF_STEP2_MUCOV,
	UKFSIF_STEP2_Q,
	UKFSIF_STEP2_chi_priory,
	UKFSIF_STEP2_x_priory,

	/*------------------------------------------------------------------------*//**
	 * @brief 	Вектор-столбец размером "Lx1"
	 */
	UKFSIF_STEP2_chi_priory_MINUS_x_priory,

	/*------------------------------------------------------------------------*//**
	 * @brief 	Вектор строка размером "1xL"
	 */
	UKFSIF_STEP2_chi_priory_MINUS_x_priory_TRANPOSE,

	/*------------------------------------------------------------------------*//**
	 * @brief 	Матрица для записи результата умножения матриц
	 * 			(chi_k|k-1 - x_k|k-1) * Transpose(chi_k|k-1 - x_k|k-1)
	 */
	UKFSIF_STEP2_RESULT_OF_MULT_2_MATRIX,

	UKFSIF_STEP2_P_apriory,

	/*------------------------------------------------------------------------*//**
	 * @brief 	Матрица размерностью "LxL*2+1"
	 * @note 	Эта матрица используется на шаге "Step2 Calculate covariance of predicted state"
	 * 			и "Step3 Calculate cross-covariance of state and output"
	 */
	UKFSIF_STEP2_chi_priory_MINUS_x_priory_TEMP,

	UKFSIF_STEP2_ARR_CELL_NUMB,
} ukfsif_step2_point_memory_place_e;

/*-------------------------------------------------------------------------*//**
 * @brief Коэффициенты для распределения сигма-точек
 */
typedef struct
{
	__UKFSIF_FPT__ alpha;
	__UKFSIF_FPT__ beta;
	__UKFSIF_FPT__ kappa;
} ukfsif_scaling_param_s;

typedef struct
{
	/*------------------------------------------------------------------------*//**
	 * @brief длина вектора пространства состояний
	 */
	uint16_t stateLen;

//	__UKFSIF_FPT__ *pMemAddr_a[UKFSIF_STEP2_MEM_CELL_NUMB];


	/*------------------------------------------------------------------------*//**
	 * @brief  Массив указателей на структуры матриц
	 */
	ukfmo_matrix_s *pMatrix_a[UKFSIF_STEP2_ARR_CELL_NUMB];

} ukfsif_step2_params_2l1_s;

//typedef struct
//{
//
//} ukfsif_step2_params_2l1_init_s;

#define ukfsif_step2_params_2l1_init_s ukfsif_step2_params_2l1_s
/*#### |End  | <-- Секция - "Определение типов" ##############################*/


/*#### |Begin| --> Секция - "Определение глобальных переменных" ##############*/
/*#### |End  | <-- Секция - "Определение глобальных переменных" ##############*/


/*#### |Begin| --> Секция - "Прототипы глобальных функций" ###################*/
extern __UKFSIF_FPT__
UKFSIF_GetLambda(
	uint16_t 		vectLen,
	__UKFSIF_FPT__ 	alpha,
	__UKFSIF_FPT__ 	kappa);

extern void
UKFSIF_InitWeightVectorMean(
	ukfsif_scaling_param_s 	*pScalParams_s,
	__UKFSIF_FPT__ 			*pWeightMean,
	uint16_t 				 vectLen);

extern void
UKFSIF_InitWeightVectorCov(
	ukfsif_scaling_param_s 	*pScalParams_s,
	__UKFSIF_FPT__ 			*pWeightCov,
	uint16_t vectLen);

extern void
UKFSIF_CalculateTheSigmaPoints_2L1(
	__UKFSIF_FPT__ *pStateVect,
	__UKFSIF_FPT__ *pSigmaPoints,
	__UKFSIF_FPT__ *pSqrtP, 			/* Указатель на двумерный массив, в котором содержится квадратный корень из матрицы ковариаций */
	__UKFSIF_FPT__  sqrtLenLambda,
	uint16_t 		stateVectLen 		/* Длина вектора пространства состояний, совпадает с количеством строк матрицы Сигма-точек */
);
/*#### |End  | <-- Секция - "Прототипы глобальных функций" ###################*/


/*#### |Begin| --> Секция - "Определение макросов" ###########################*/
/*#### |End  | <-- Секция - "Определение макросов" ###########################*/


/*#### |Begin| --> Секция - "Include - подмодули" ############################*/
/*#### |End  | <-- Секция - "Include - подмодули" ############################*/

#endif	/* LIB_A_UKFSIF_UKF_STANDART_INIT_FNC_H_ */

/*############################################################################*/
/*################################ END OF FILE ###############################*/
/*############################################################################*/
