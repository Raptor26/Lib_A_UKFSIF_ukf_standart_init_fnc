/**
 * @file   	%<%NAME%>%.%<%EXTENSION%>%
 * @author 	%<%USER%>%
 * @version
 * @date 	%<%DATE%>%, %<%TIME%>%
 * @brief 	Правила присваивания имен матрицам:
 *			"k-1"	- "predict"
 *			"k|k-1"	- "apriori"
 *			"k"		- "posteriori"
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
	UKFSIF_INIT_muMean = 0u,
	UKFSIF_INIT_muCovar,
	UKFSIF_INIT_Q,
	UKFSIF_INIT_R,
	UKFSIF_INIT_P,
	UKFSIF_INIT_P_SQRT,
	UKFSIF_INIT_P_apriory,

	UKFSIF_INIT_chi_predict,
	UKFSIF_INIT_chi_apriori,

	UKFSIF_INIT_x_apriori,

	UKFSIF_INIT_STEP2_chi_priory_MINUS_x_priory,

	UKFSIF_INIT_STEP2_chi_priory_MINUS_x_priory_TRANSPOSE,

	/* @todo получится ли использовать эту матрицу на следующих шагах:
	 * - Calculate covariance of predicted state
	 * - Calculate covariance of predicted output
	 * - Calculate cross-covariance of state and output */
	UKFSIF_INIT_STEP2_result_of_mult_2_matrix,

	UKFSIF_INIT_ARR_CELL_NUMB,
} ukfsif_init_point_mem_matrix_struct_e;

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

typedef enum
{
	UKFSIF_STEP3_MUMEAN = 0,
	UKFSIF_STEP3_MUCOV,

	UKFSIF_STEP3_R,

	UKFSIF_STEP3_psi_apriori,

	UKFSIF_STEP3_y_apriori,


	UKFSIF_STEP3_ARR_CELL_NUMB,
} ukfsif_step3_point_memory_place_e;

typedef enum
{
	/*------------------------------------------------------------------------*//**
	 * @brief Вектор-столбец для хранения весовых коэффициентов
	 *
	 * @note  	Размерность: "(2L+1)xL"
	 */
	UKFSIF_CALC_COVAR_GENERIC_muCovar = 0,

	/*------------------------------------------------------------------------*//**
	 * @brief 	Матрица преобразованных Сигма-точек:
	 *         	- "chi_k|k-1"
	 *         	- "psi_k|k-1"
	 *
	 * @note  	Размерность: "Lx(2L+1)"
	 */
	UKFSIF_CALC_COVAR_GENERIC_sigma_apriori,

	/*------------------------------------------------------------------------*//**
	 * @brief  	Матрица усредненной матрицы Сигма-точек
	 *          - "x_k|k-1"
	 *          - "y_k|k-1"
	 *
	 * @note  	Размерность: "Lx1"
	 */
	UKFSIF_CALC_COVAR_GENERIC_vect_apriori,

	/*------------------------------------------------------------------------*//**
	 * @brief  Вектор-столбец для хранения:
	 *         - "chi_k|k-1" - "x_k|k-1"
	 *         - "psi_k|k-1" - "y_k|k-1"
	 *
	 * @note  	Размерность: "Lx1"
	 */
	UKFSIF_CALC_COVAR_GENERIC_sigma_apriori_MINUS_state_apriori,

	/*------------------------------------------------------------------------*//**
	 * @brief  Вектор-столбец для хранения:
	 *         - Transpose("chi_k|k-1" - "x_k|k-1")
	 *         - Transpose("psi_k|k-1" - "y_k|k-1")
	 *
	 * @note  	Размерность: "Lx1"
	 */
	UKFSIF_CALC_COVAR_GENERIC_sigma_apriori_MINUS_state_apriori_TRANSPOSE,

	/*------------------------------------------------------------------------*//**
	 * @brief  Матрица для хранения результата умножения:
	 *         - ("chi_k|k-1" - "x_k|k-1") * Transpose("chi_k|k-1" - "x_k|k-1")
	 *         - ("psi_k|k-1" - "y_k|k-1") * Transpose("psi_k|k-1" - "y_k|k-1")
	 *
	 * @note 	Размерность: "LxL"
	 */
	UKFSIF_CALC_COVAR_GENERIC_MULT_2_MATRIX,

	/*------------------------------------------------------------------------*//**
	 * @brief 	Матрица ковариации, сюда будет записан результат
	 *
	 * @note  	Размерность: "LxL"
	 */
	UKFSIF_CALC_COVAR_GENERIC_covariance_apriori,

	/*------------------------------------------------------------------------*//**
	 * @brief  Размер массива указателей на структуры данных
	 */
	UKFSIF_CALC_COVAR_GENERIC_ARR_CELL_NUMB,
} ukfsif_calc_covar_generic_e;

typedef enum
{
	/*------------------------------------------------------------------------*//**
	 * @brief  Вектор-столбец для хранения:
	 *         - "chi_k|k-1" - "x_k|k-1"
	 *         - "psi_k|k-1" - "y_k|k-1"
	 *
	 * @note  	Размерность: "Lx(2L+1)"
	 */
	UKFSIF_CALC_MEAN_GENERIC_sigma_apriori = 1u,

	/*------------------------------------------------------------------------*//**
	 * @brief Вектор-столбец для хранения весовых коэффициентов
	 *
	 * @note  	Размерность: "(2L+1)x1"
	 */
	UKFSIF_CALC_MEAN_GENERIC_muMean,

	/*------------------------------------------------------------------------*//**
	 * @brief Вектор-столбец для записи усредненного значения Сигма-точек
	 *
	 * @note  Размерность: "Lx1"
	 */
	UKFSIF_CALC_MEAN_GENERIC_vect_apriori,

	/*------------------------------------------------------------------------*//**
	 * @brief  Размер массива указателей на структуры данных
	 */
	UKFSIF_CALC_MEAN_GENERIC_ARR_CELL_NUMB,
} ukfsif_calc_mean_generic_e;

typedef enum
{
	UKFSIF_CALC_KALMAN_GAIN_Pxy = 0u,
	UKFSIF_CALC_KALMAN_GAIN_Pyy,
	UKFSIF_CALC_KALMAN_GAIN_Pyy_INV,
	UKFSIF_CALC_KALMAN_GAIN_K,

	UKFSIF_CALC_KALMAN_GAIN_ARR_CELL_NUMB,
} ukfsif_calc_kalman_gain_e;

typedef enum
{
	/*------------------------------------------------------------------------*//**
	 * @brief Вектор-столбец пространства состояний после проведения коррекции
	 * 
	 * @note  Размерность: "Lx1"
	 */
	UKFSIF_UPDATE_STATE_ESTIMATE_x_posteriori = 0u,

	/*------------------------------------------------------------------------*//**
	 * @brief Вектор-столбец пространства состояний до проведения коррекции
	 * 
	 * @note  Размерность: "Lx1"
	 */
	UKFSIF_UPDATE_STATE_ESTIMATE_x_apriori,

	/*------------------------------------------------------------------------*//**
	 * @brief  Матрица коэффициентов усиления Калмана
	 * 
	 * @note  Размерность: "LxL"
	 */
	UKFSIF_UPDATE_STATE_ESTIMATE_K,

	/*------------------------------------------------------------------------*//**
	 * @brief  Вектор-столбец оценки вектор-столбца измерения 
	 * 
	 * @note  Размерность: "Lx1"
	 */
	UKFSIF_UPDATE_STATE_ESTIMATE_y_apriori,

	/*------------------------------------------------------------------------*//**
	 * @brief  Вектор-столбец измерений
	 * 
	 * @note  Размерность: "Lx1"
	 */
	UKFSIF_UPDATE_STATE_ESTIMATE_meas,

	UKFSIF_UPDATE_STATE_ESTIMATE_innovation,

	/*------------------------------------------------------------------------*//**
	 * @brief  Размер массива указателей на структуры данных
	 */
	UKFSIF_UPDATE_STATE_ESTIMATE_ARR_CELL_NUMB,
} ukfsif_update_state_estimate_e;

typedef enum
{
	/*------------------------------------------------------------------------*//**
	 * @brief  
	 * 
	 * @note  Размерность: "LxL"
	 */
	UKFSIF_UPDATE_ERR_COVAR_P_posteriori = 0u,

	/*------------------------------------------------------------------------*//**
	 * @brief  
	 * 
	 * @note  Размерность: "LxL"
	 */
	UKFSIF_UPDATE_ERR_COVAR_P_apriori,

	/*------------------------------------------------------------------------*//**
	 * @brief  Матрица коэффициентов усиления 
	 * 
	 * @note  Размерность: "LxL"
	 */
	UKFSIF_UPDATE_ERR_COVAR_K,

		/*------------------------------------------------------------------------*//**
	 * @brief  
	 * 
	 * @note  Размерность: "LxL"
	 */
	UKFSIF_UPDATE_ERR_COVAR_K_TRANSPOSE,

	/*------------------------------------------------------------------------*//**
	 * @brief  
	 * 
	 * @note  Размерность: "LxL"
	 */
	UKFSIF_UPDATE_ERR_COVAR_Pyy,

	UKFSIF_UPDATE_ERR_COVAR_ARR_CELL_NUMB,
} ukfsif_update_err_covar_e;

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

	/*------------------------------------------------------------------------*//**
	 * @brief  Массив указателей на структуры матриц
	 */
	ukfmo_matrix_s *pMatrix_a[UKFSIF_STEP2_ARR_CELL_NUMB];

} ukfsif_step2_params_2l1_s;

typedef struct
{
	/*------------------------------------------------------------------------*//**
	 * @brief длина вектора пространства состояний
	 */
	uint16_t stateLen;

	ukfmo_matrix_s *pMatrix_a[UKFSIF_STEP3_ARR_CELL_NUMB];
} ukfsif_step3_params_2l1_s;

/*-------------------------------------------------------------------------*//**
 * @brief  Структура для хранения указателей на матрицы, необходимые для
 *         расчета ковариации
 */
typedef struct
{
	/*------------------------------------------------------------------------*//**
	 * @brief Длина вектора пространства состояний
	 */
	uint16_t stateLen;

	/*------------------------------------------------------------------------*//**
	 * @brief  Массив указателей на структуры, содержащие матрицы
	 */
	ukfmo_matrix_s *pMatrix_a[UKFSIF_CALC_COVAR_GENERIC_ARR_CELL_NUMB];
} ukfsif_calc_covar_generic_s;

/*-------------------------------------------------------------------------*//**
 * @brief  Структура для расчета "среднего" от матрицы Сигма-точек
 * 
 *         (для преобразования матрицы Сигма-точек размерностью Lx(2L+1))
 *         в вектор-столбец размерностью Lx1 с помощью вектора весовых 
 *         коэффициентов размерностью (2L+1)x1)
 */
typedef struct
{
	uint16_t stateLen;

	ukfmo_matrix_s *pMatrix_a[UKFSIF_CALC_MEAN_GENERIC_ARR_CELL_NUMB];
} ukfsif_calc_mean_generic_s;

typedef struct
{
	uint16_t stateLen;

	ukfmo_matrix_s *pMatrix_a[UKFSIF_CALC_KALMAN_GAIN_ARR_CELL_NUMB];
} ukfsif_calc_kalman_gain_s;

typedef struct
{
	uint16_t stateLen;

	ukfmo_matrix_s *pMatrix_a[UKFSIF_UPDATE_STATE_ESTIMATE_ARR_CELL_NUMB];
} ukfsif_update_state_s;

typedef struct
{
	uint16_t stateLen;

	ukfmo_matrix_s *pMatrix_a[];
} ukfsif_update_err_covar_s;

typedef struct
{
	ukfsif_update_err_covar_s 	updateErrCov_s;

	ukfsif_update_state_s 		updateState_s;

	ukfsif_calc_kalman_gain_s 	calcKalmanGain_s;
} ukfsif_all_data_s;

typedef struct
{

} ukfsif_all_data_init_s;

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
