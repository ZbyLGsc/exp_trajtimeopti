#include <pcl/io/pcd_io.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/registration/gicp.h>
#include <pcl/visualization/cloud_viewer.h>
#include "gicp/gicp.h"


pcl::PointCloud<pcl::PointXYZRGBA>::Ptr modelCloud(new pcl::PointCloud<pcl::PointXYZRGBA>);
pcl::PointCloud<pcl::PointXYZRGBA>::Ptr modelCloudDownsampled(new pcl::PointCloud<pcl::PointXYZRGBA>);
pcl::PointCloud<pcl::PointXYZRGBA>::Ptr dataCloud(new pcl::PointCloud<pcl::PointXYZRGBA>);
pcl::PointCloud<pcl::PointXYZRGBA>::Ptr dataCloudDownsampled(new pcl::PointCloud<pcl::PointXYZRGBA>);

pcl::PointCloud<pcl::PointXYZRGBA>::Ptr transformed(new pcl::PointCloud<pcl::PointXYZRGBA>);

void viewerOneOff(pcl::visualization::PCLVisualizer& viewer)
{
    viewer.addPointCloud(modelCloudDownsampled, "model");
    //viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 255, 0, 0, "model");
    viewer.addPointCloud(transformed, "transformed");
    //viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 0, 255, "transformed");
}


int main(int argv, char **args)
{
    // load files and convert to point type pcl::PointXYZRGBA
    sensor_msgs::PointCloud2 cloudBlob;

    pcl::io::loadPCDFile("1311868164.399026009.pcd", cloudBlob);
    pcl::fromROSMsg(cloudBlob, *modelCloud);
    pcl::io::loadPCDFile("1311868164.698987188.pcd", cloudBlob);
    pcl::fromROSMsg(cloudBlob, *dataCloud);

    // make dense
    std::vector<int> indices;
    pcl::removeNaNFromPointCloud(*modelCloud, *modelCloud, indices);
    pcl::removeNaNFromPointCloud(*dataCloud,  *dataCloud,  indices);

    // downsample clouds
    pcl::VoxelGrid<pcl::PointXYZRGBA> vg;
    vg.setInputCloud(modelCloud);
    vg.setLeafSize (0.01f, 0.01f, 0.01f);
    vg.filter(*modelCloudDownsampled);

    vg.setInputCloud(dataCloud);
    vg.setLeafSize (0.01f, 0.01f, 0.01f);
    vg.filter(*dataCloudDownsampled);


    /* Segal-GICP */
    dgc::gicp::GICPPointSet ps_in, ps_out;
    
    for(int i=0; i<dataCloudDownsampled->points.size(); ++i)
    {
        dgc::gicp::GICPPoint p;
        p.x = dataCloudDownsampled->points[i].x;
        p.y = dataCloudDownsampled->points[i].y;
        p.z = dataCloudDownsampled->points[i].z;
        
        ps_in.AppendPoint(p);
    }
    
    for(int i=0; i<modelCloudDownsampled->points.size(); ++i)
    {
        dgc::gicp::GICPPoint p;
        p.x = modelCloudDownsampled->points[i].x;
        p.y = modelCloudDownsampled->points[i].y;
        p.z = modelCloudDownsampled->points[i].z;
        
        ps_out.AppendPoint(p);
    }    
    
    dgc_transform_t T_guess, T_delta, T_est;
    dgc_transform_identity(T_delta);
    dgc_transform_identity(T_guess);

    ps_in.BuildKDTree();
    ps_out.BuildKDTree();
    ps_in.ComputeMatrices();
    ps_out.ComputeMatrices();
    ps_out.SetDebug(0);
    //ps_out.SetMaxIterationInner(8);
   
    ps_out.AlignScan(&ps_in, T_guess, T_delta, 5);
    dgc_transform_copy(T_est, T_guess);
    dgc_transform_left_multiply(T_est, T_delta);          
    
    Eigen::Matrix4f em;
    
    for(int i=0; i<4; ++i)
        for(int j=0; j<4; ++j)
            em(i,j) = (float) T_est[i][j];
            
    std::cout << em << std::endl << std::endl;


    // PCL-GICP
    pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZRGBA, pcl::PointXYZRGBA> gicp;
    gicp.setInputCloud(dataCloudDownsampled);
    gicp.setInputTarget(modelCloudDownsampled);

    gicp.align(*transformed);    

    //std::cout << "has converged: " << gicp.hasConverged() << " score: " << gicp.getFitnessScore() << std::endl << std::endl;
    std::cout << gicp.getFinalTransformation() << std::endl;
    
    // apply transform to data cloud (to fit model cloud)
    pcl::transformPointCloud(*dataCloudDownsampled, *transformed, gicp.getFinalTransformation());

    // show results
    pcl::visualization::CloudViewer viewer("Simple Cloud Viewer");
    viewer.runOnVisualizationThreadOnce(viewerOneOff);

    while(!viewer.wasStopped())
    {
    }

    return 0;
}
