/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DataSet_h
#define _DataSet_h

#include "itkObject.h"
#include "itkObjectFactory.h"

#include <iostream>
#include <string>


class DataSet : public itk::Object
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	typedef DataSet Self;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	itkNewMacro(Self);

	// Make sure all values are OK
	virtual bool CheckValues();

    virtual void PrintSelf(std::ostream& os);

	// Name of the XML file that has been used to generate this DataSet object
	inline void SetXMLFileName(const char* fn){ m_XMLFileName = fn; }
	inline const char* GetXMLFileName() const { return m_XMLFileName; }
	
	
	// Compulsory parameters
	itkGetMacro(DeformableObjectType, std::string);
	itkSetMacro(DeformableObjectType, std::string);

	itkGetMacro(DataSigma, double);
	itkSetMacro(DataSigma, double);


protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	DataSet();

	~DataSet();


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	// Name of the XML file that has been used to generate this DataSet object
	const char * m_XMLFileName;
	
	
	std::string m_DeformableObjectType;

	double m_DataSigma;
	/******************************************************************************************************************************
	                                                               CHANGE BEGIN
	******************************************************************************************************************************/
	double m_DataSigma_Normalized_Hyperparameter;
	double m_DataSigma_Prior;
	/******************************************************************************************************************************
	                                                           CHANGE END
	******************************************************************************************************************************/

	std::string m_KernelType;

	double m_KernelWidth;
	// double m_P3MWorkingSpacingRatio;
	// double m_P3MPaddingFactor;

	double m_ImageGridDownsampling;

	std::string m_AnatomicalCoordinateSystem;

	bool m_reOrient;


}; /* class DataSet */


#endif /* _DataSet_h */
