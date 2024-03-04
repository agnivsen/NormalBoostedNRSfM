# Normal Boosted Non-Rigid Structure-from-Motion
Code and dataset associated with the [IEEE ICRA 2024](https://2024.ieee-icra.org/) publication 'Using Specularities to Boost Non-Rigid Structure-from-Motion'

In this article and this GitHub repository, we present a method to reconstruct the time-varying shape of a non-rigid/deforming surface from the following inputs:

   1. Accurate keypoint correspondences (e.g.: from feature matching from standard feature descriptors) on the deforming surface captured across multiple images (preferably more than 5 images)

   2. Sparse surface normals on the surface of the object (**optional input to the code**)


Despite promising use-cases such as the grasping of deformable objects and visual navigation in a non-rigid environment, Non-Rigid Structure-from-Motion (NRSfM) has had limited applications in robotics due to a lack of accuracy. To remedy this, we propose a new method which boosts the accuracy of NRSfM using sparse surface normals. Surface normal information is available from many sources, including structured lighting, homography decomposition of infinitesimal planes and shape priors. However, these sources are not always available. We thus propose a widely available new source of surface normals: the specularities. Our first technical contribution is a method which detects specular highlights and reconstructs the surface normals from it. It assumes that the light source is approximately localised, which is widely applicable in robotics applications such as endoscopy. Our second technical contribution is an NRSfM method which exploits a sparse surface normal set. For that, we propose a novel convex formulation and a globally optimal solution method.

A descriptive YouTube video is available here:

[<img src="https://github.com/agnivsen/NormalBoostedNRSfM/assets/5153445/7d0a725d-62b2-462a-a51a-289ab0c910d9">](https://www.youtube.com/watch?v=jNhC3noMyEs)

More details about the paper along with a copy of this code and additional dataset are available at the website of the EnCoV research group: https://encov.ip.uca.fr/ab/code_and_datasets/

# Dependencies

This code depends on the following external libraries/toolboxes:

 -  [CVX](http://cvxr.com/cvx): the MATLAB software for disciplined convex programming
 -  Developable surface simulator from [Perriollat et al., 2013]

# Data format

All data have been provided in standard NRSfM/SfT format, we explain it below:

Say the data contains _n_ images tracking up to _m_ feature correspondences across images.

The data is a MATLAB _struct_, say **Data**. It contains the following subfields:

* Data.Pgth(i).P is a [3 x _m_] matrix, containing 3D groundtruth points, for all _i_ in [1, _n_]

* Data.p(i).p is a [2 x _m_] or [3 x _m_] matrix, containing tracked point correspondences (image coordinates), for all _i_ in [1, _n_]. If the matrix is of size [3 x _m_], the last row is necessarily all ones

* Data.v is a [_n_ x _m_] matrix which is zero if a point (indexed by column) is invisible at that particular image corresponding to the row number, one otherwise.



# Scripts

Run the following script:

 - **runMe.m:** invokes our proposed NRSfM method with either synthetic data generated using [Perriollat et al., 2013] or the real data introduced in our paper

# Citation

This article has been accepted at the IEEE International Conference of Robotics and Automation (ICRA) 2024. If you find this code and the associated paper useful, you may cite us with:
```
@inproceedings{sengupta2024specular,
  title={Using Specularities to Boost Non-Rigid Structure-from-Motion},
  author={Sengupta, Agniva and Makki, Karim and Bartoli, Adrien},
  booktitle={2024 IEEE International Conference on Robotics and Automation (ICRA)},
  year={2024},
  organization={IEEE}
}
```
|[EnCoV (PrePrint)](http://encov.ip.uca.fr/publications/pubfiles/2024_Sengupta_etal_ICRA_normal.pdf)|


---
## References


[Perriollat et al., 2013]: Perriollat, M., & Bartoli, A. (2013). A computational model of bounded developable surfaces with application to image‐based three‐dimensional reconstruction. Computer Animation and Virtual Worlds, 24(5), 459-476.


