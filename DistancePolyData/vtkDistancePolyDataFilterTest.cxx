#include <vtkVersion.h>
#include <vtkSmartPointer.h>
 
#include <vtkActor.h>
#include <vtkDistancePolyDataFilter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkCleanPolyData.h>
#include <vtkPointData.h>

#include <cmath>
#include <cstdlib>

void printUsage(int argc, char** argv)
{
  std::cerr << "Call: " << argc << std::endl;
  std::cerr << argv[0] << " args : inputMesh1 inputMesh2 outputMesh1 outputMesh2 "
                       << "[weightsForDistanceMesh] [weightsForDistanceThreshold=0.1]"
             << std::endl;
  std::cerr << "Takes two meshes "
            << "and outputs mesh with point-to-cell distances as scalars."
            << std::endl
            << "If weightsForDistanceMesh specified, the mean distance over voxels "
            << "where the weightsForDistanceMesh scalars greater than weightsForDistanceThreshold."
            << std::endl;
}

int main(int argc, char* argv[])
{
    ////////////
    // Command Line Arguments
    if (argc < 5)
    {
        printUsage(argc, argv);
        return EXIT_FAILURE;
    }

    vtkSmartPointer<vtkPolyData> input1;
    vtkSmartPointer<vtkPolyData> input2;
    vtkSmartPointer<vtkPolyData> input3;
    vtkSmartPointer<vtkPolyDataReader> reader1 =
      vtkSmartPointer<vtkPolyDataReader>::New();
    reader1->SetFileName(argv[1]);
    reader1->Update();
    input1 = reader1->GetOutput();
 
    vtkSmartPointer<vtkPolyDataReader> reader2 =
      vtkSmartPointer<vtkPolyDataReader>::New();
    reader2->SetFileName(argv[2]);
    reader2->Update();
    input2 = reader2->GetOutput();

    vtkSmartPointer<vtkPolyDataReader> reader3 =
      vtkSmartPointer<vtkPolyDataReader>::New();
    if(argc>5)
    {
      reader3->SetFileName(argv[5]);
      reader3->Update();
      input3 = reader3->GetOutput();
    }

    double weightsForMeanThreshold = 0.1;
    if(argc>6)
    {
      weightsForMeanThreshold = std::atof(argv[6]);
    }

    /// check if same number of points
    vtkSmartPointer<vtkPoints> shapePoints1 = input1->GetPoints();
    vtkIdType shapePointNo1 = shapePoints1->GetNumberOfPoints();

    if(argc>5)
    {
      vtkIdType shapePointNo3 = input3->GetPoints()->GetNumberOfPoints();

      if(shapePointNo1 != shapePointNo3)
      {
        std::cerr << "ERROR: The inputMesh1 do not have the same "
                  << "number of points as weightsForMeanMesh!" << std::endl;
        return 1;
      }
    }


    vtkSmartPointer<vtkCleanPolyData> clean1 =
      vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
    clean1->SetInputConnection( input1->GetProducerPort());
#else
    clean1->SetInputData( input1);
#endif
 
    vtkSmartPointer<vtkCleanPolyData> clean2 =
      vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
    clean2->SetInputConnection( input2->GetProducerPort());
#else
    clean2->SetInputData( input2);
#endif
 
    vtkSmartPointer<vtkDistancePolyDataFilter> distanceFilter =
      vtkSmartPointer<vtkDistancePolyDataFilter>::New();
 
    distanceFilter->SetInputConnection( 0, clean1->GetOutputPort() );
    distanceFilter->SetInputConnection( 1, clean2->GetOutputPort() );
    distanceFilter->NegateDistanceOn();
    distanceFilter->Update();
 
    /*
     * Save the first output
     */
    vtkSmartPointer<vtkPolyDataWriter> writer
      = vtkSmartPointer<vtkPolyDataWriter>::New();
#if VTK_MAJOR_VERSION <= 5
    writer->SetInputConnection(distanceFilter->GetOutput());
#else
    writer->SetInputData(distanceFilter->GetOutput());
#endif
    writer->SetFileName(argv[3]);
    writer->Update();

    /// iterate through points
    double meanDistance = 0.0;
    unsigned long count = 0;

    for (vtkIdType ii = 0; ii < shapePointNo1; ++ii)
    {
      if(argc<=5 ||
              input3->GetPointData()->GetScalars()->GetTuple1(ii)>=weightsForMeanThreshold)
      {
        meanDistance += std::abs(distanceFilter->GetOutput()->GetPointData()->GetScalars()->GetTuple1(ii));
        ++count;
      }
    }

    if(count)
    {
      std::cout << "Mean distance over specified region: " << meanDistance/count
                << " (" << count << " points)." << std::endl;
    }
    else
    {
      std::cout << "No points found in the specified region for mean distance"
                  << std::endl;
    }




    // Compute the distance in other direction
    vtkSmartPointer<vtkDistancePolyDataFilter> distanceFilter2 =
      vtkSmartPointer<vtkDistancePolyDataFilter>::New();

    distanceFilter2->SetInputConnection( 1, clean1->GetOutputPort() );
    distanceFilter2->SetInputConnection( 0, clean2->GetOutputPort() );
    distanceFilter2->Update();

    /*
     * Save the second output
     */
    vtkSmartPointer<vtkPolyDataWriter> writer2
      = vtkSmartPointer<vtkPolyDataWriter>::New();
#if VTK_MAJOR_VERSION <= 5
    writer2->SetInputConnection(distanceFilter2->GetOutput());
#else
    writer2->SetInputData(distanceFilter2->GetOutput());
#endif
    writer2->SetFileName(argv[4]);
    writer2->Update();

    /// iterate through points
    double meanDistance2 = 0.0;
    unsigned long count2 = 0;

    for (vtkIdType ii = 0; ii < shapePointNo1; ++ii)
    {
      if(argc<=5 ||
              input3->GetPointData()->GetScalars()->GetTuple1(ii)>=weightsForMeanThreshold)
      {
        meanDistance2 += std::abs(distanceFilter2->GetOutput()->GetPointData()->GetScalars()->GetTuple1(ii));
        ++count2;
      }
    }

    if(count2)
    {
      std::cout << "Mean distance over specified region: " << meanDistance2/count2
                << " (" << count2 << " points)." << std::endl;
    }
    else
    {
      std::cout << "No points found in the specified region for mean distance"
                  << std::endl;
    }

    return EXIT_SUCCESS;
}
