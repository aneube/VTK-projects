#include <vtkVersion.h>
#include <vtkSmartPointer.h>
 
#include <vtkActor.h>
#include <vtkDistancePolyDataFilter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkCleanPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>

void printUsage(int argc, char** argv)
{
  std::cerr << "Call: " << argc << std::endl;
  std::cerr << argv[0] << " args : inputMesh1 inputMesh2 outputMesh1 outputMesh2 "
                       << "[weightsForDistanceMesh1] [weightsForDistanceThreshold1=0.1]"
                       << "[weightsForDistanceMesh2] [weightsForDistanceThreshold2=0.1]"
             << std::endl;
  std::cerr << "Takes two meshes "
            << "and outputs mesh with point-to-cell distances as scalars."
            << std::endl
            << "If weightsForDistanceMesh(1,2) specified, the mean distance over voxels "
            << "where the weightsForDistanceMesh(1,2) scalars greater than weightsForDistanceThreshold(1,2)."
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
    vtkSmartPointer<vtkPolyData> input4;
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

    vtkSmartPointer<vtkPolyDataReader> reader4 =
      vtkSmartPointer<vtkPolyDataReader>::New();
    if(argc>7)
    {
      reader4->SetFileName(argv[7]);
      reader4->Update();
      input4 = reader4->GetOutput();
    }
    double weightsForMeanThreshold2 = 0.1;
    if(argc>8)
    {
      weightsForMeanThreshold2 = std::atof(argv[8]);
    }


    // Prepare the text file to save all the distances
    std::string output1_filename = argv[3];
    std::string distances_filename = output1_filename.replace(output1_filename.end()-4, output1_filename.end(), "_distances.txt");
    ofstream distances_file;
    distances_file.open(distances_filename.c_str());

    /// check if same number of points
    vtkSmartPointer<vtkPoints> shapePoints1 = input1->GetPoints();
    vtkIdType shapePointNo1 = shapePoints1->GetNumberOfPoints();
    vtkSmartPointer<vtkPoints> shapePoints2 = input2->GetPoints();
    vtkIdType shapePointNo2 = shapePoints2->GetNumberOfPoints();

    if(argc>5)
    {
      vtkIdType shapePointNo3 = input3->GetPoints()->GetNumberOfPoints();

      if(shapePointNo1 != shapePointNo3)
      {
        std::cerr << "ERROR: The inputMesh1 do not have the same "
                  << "number of points as weightsForMeanMesh1!" << std::endl;
        return 1;
      }
    }

    if(argc>7)
    {
      vtkIdType shapePointNo4 = input4->GetPoints()->GetNumberOfPoints();

      if(shapePointNo2 != shapePointNo4)
      {
        std::cerr << "ERROR: The inputMesh2 do not have the same "
                  << "number of points as weightsForMeanMesh2!" << std::endl;
        return 1;
      }
    }

/*
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
*/
    vtkSmartPointer<vtkDistancePolyDataFilter> distanceFilter =
      vtkSmartPointer<vtkDistancePolyDataFilter>::New();
 
    //distanceFilter->SetInputConnection( 0, clean1->GetOutputPort() );
    //distanceFilter->SetInputConnection( 1, clean2->GetOutputPort() );
    distanceFilter->SetInputData( 0, input1 );
    distanceFilter->SetInputData( 1, input2 );

    distanceFilter->NegateDistanceOn();
    distanceFilter->Update();
 
    vtkSmartPointer<vtkPolyData> output1 = distanceFilter->GetOutput();

    /*
     * Compute stats
     */
    /// iterate through points
    vtkSmartPointer<vtkDoubleArray> scalars =
        vtkSmartPointer<vtkDoubleArray>::New();
     scalars->SetNumberOfValues(shapePointNo1);

    double meanDistance = 0.0;
    double meanDistanceSquared = 0.0;
    unsigned long count = 0;
    double distance = 0.0;
    std::vector<double> distance_values;

    for (vtkIdType ii = 0; ii < shapePointNo1; ++ii)
    {
      if(argc<=5 ||
              input3->GetPointData()->GetScalars()->GetTuple1(ii)>=weightsForMeanThreshold)
      {
        scalars->SetValue(ii, output1->GetPointData()->GetScalars()->GetTuple1(ii));
        distance = std::abs(output1->GetPointData()->GetScalars()->GetTuple1(ii));
        distance_values.push_back(distance);
        meanDistance += distance;
        meanDistanceSquared += distance * distance;
        ++count;

        distances_file << distance << "\n";
      }
      else
      {
        scalars->SetValue(ii, 0);
      }
    }

    double sigma;
    double conf_interval;
    if(count)
    {
      // Compute the standard deviation and the 95% confidence interval
      meanDistance /= count;
      meanDistanceSquared /= count;
      sigma = std::sqrt(meanDistanceSquared - meanDistance * meanDistance);
      conf_interval = 1.96 * sigma / std::sqrt(static_cast<double>(count));

      std::cout << "Mean distance over specified region: " << meanDistance
                << " +- " << sigma << " (conf.int +- " << conf_interval
                << " " << count << " points)." << std::endl;
    }
    else
    {
      std::cout << "No points found in the specified region for mean distance"
                  << std::endl;
    }
    output1->GetPointData()->SetScalars(scalars);

    /*
     * Save the first output
     */
    vtkSmartPointer<vtkPolyDataWriter> writer
      = vtkSmartPointer<vtkPolyDataWriter>::New();
#if VTK_MAJOR_VERSION <= 5
    writer->SetInputConnection(output1);
#else
    writer->SetInputData(output1);
#endif
    writer->SetFileName(argv[3]);
    writer->Update();



    // Compute the distance in other direction
    vtkSmartPointer<vtkDistancePolyDataFilter> distanceFilter2 =
      vtkSmartPointer<vtkDistancePolyDataFilter>::New();

//    distanceFilter2->SetInputConnection( 1, clean1->GetOutputPort() );
//    distanceFilter2->SetInputConnection( 0, clean2->GetOutputPort() );
    distanceFilter2->SetInputData( 0, input2 );
    distanceFilter2->SetInputData( 1, input1 );


    distanceFilter2->Update();

    vtkSmartPointer<vtkPolyData> output2 = distanceFilter2->GetOutput();

    /*
     * Compute stats
     */
    /// iterate through points
    vtkSmartPointer<vtkDoubleArray> scalars2 =
        vtkSmartPointer<vtkDoubleArray>::New();
    scalars2->SetNumberOfValues(shapePointNo2);

    /// iterate through points
    double meanDistance2 = 0.0;
    double meanDistance2Squared = 0.0;
    unsigned long count2 = 0;

    for (vtkIdType ii = 0; ii < shapePointNo2; ++ii)
    {
      if(argc<=7 ||
              input4->GetPointData()->GetScalars()->GetTuple1(ii)>=weightsForMeanThreshold2)
      {
        scalars2->SetValue(ii, output2->GetPointData()->GetScalars()->GetTuple1(ii));
        distance = std::abs(output2->GetPointData()->GetScalars()->GetTuple1(ii));
        distance_values.push_back(distance);
        meanDistance2 += distance;
        meanDistance2Squared += distance * distance;
        ++count2;

        distances_file << distance << "\n";
      }
      else
      {
        scalars2->SetValue(ii, 0);
      }
    }

    double sigma2;
    double conf_interval2;
    if(count2)
    {
      // Compute the standard deviation and the 95% confidence interval
      meanDistance2 /= count2;
      meanDistance2Squared /= count2;
      sigma2 = std::sqrt(meanDistance2Squared - meanDistance2 * meanDistance2);
      conf_interval2 = 1.96 * sigma2 / std::sqrt(static_cast<double>(count2));

      std::cout << "Mean distance over specified region: " << meanDistance2
                << " +- " << sigma2 << " (conf.int +- " << conf_interval2
                << " " << count2 << " points)." << std::endl;
    }
    else
    {
      std::cout << "No points found in the specified region for mean distance"
                  << std::endl;
    }
    output2->GetPointData()->SetScalars(scalars2);

    if(count && count2)
    {
        // compute quantiles
        std::sort(distance_values.begin(), distance_values.end());
        int size = distance_values.size();
        int q90 = static_cast<int>(0.9*size + 0.5);
        int q95 = static_cast<int>(0.95*size + 0.5);

        std::cout << "Average mean distance over specified region: " << 0.5*(meanDistance + meanDistance2)
                  << " +- " << 0.5 * (sigma + sigma2) << " , conf.int +- " << 0.5 * (conf_interval + conf_interval2)
                  << " : " << 0.5*(meanDistance + meanDistance2) - 0.5 * (conf_interval + conf_interval2)
                  << " - " << 0.5*(meanDistance + meanDistance2) + 0.5 * (conf_interval + conf_interval2)
                  << " , 90% quantile q90: " << distance_values[q90]
                  << " , 95% quantile q95: " << distance_values[q95]
                  << " , max: " << distance_values[size-1]
                  << std::endl;
    }

    /*
     * Save the second output
     */
    vtkSmartPointer<vtkPolyDataWriter> writer2
      = vtkSmartPointer<vtkPolyDataWriter>::New();
#if VTK_MAJOR_VERSION <= 5
    writer2->SetInputConnection(output2);
#else
    writer2->SetInputData(output2);
#endif
    writer2->SetFileName(argv[4]);
    writer2->Update();

    distances_file.close();

    return EXIT_SUCCESS;
}
