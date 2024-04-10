from typing import Literal
import numpy as np
import numpy.typing as npt
import random


def squaredErrorDistortion(data: npt.NDArray, centers: npt.NDArray):
    n, m = data.shape
    k = centers.shape[0]
    min_dist_to_centers = minDistanceToCenters(data, centers)

    distortion = 1/n * np.sum(min_dist_to_centers * min_dist_to_centers)
    return distortion


def distanceToCenters(data: npt.NDArray, centers: npt.NDArray):
    n, m = data.shape
    k = centers.shape[0]
    distances = np.zeros((n, k), dtype=data.dtype)
    for i, center in enumerate(centers):
        dist = np.linalg.norm(data-center, axis=1)
        distances[:, i] = dist
    return distances

def argminDistanceToCenters(data: npt.NDArray, centers: npt.NDArray):
    """Returns an array with the minium distance of each point in data to all centers

    Args:
        data (npt.NDArray): n x m array
        centers (npt.NDArray): k x m array

    Returns:
        argmin_dist_to_centers (npt.NDArray): 1-D array of lenght n
    """
    distances = distanceToCenters(data, centers)
    closest_centers = np.argmin(distances, axis=1)
    return closest_centers


def minDistanceToCenters(data: npt.NDArray, centers: npt.NDArray):
    """Returns an array with the minimum distance of each point in data to all centers

    Args:
        data (npt.NDArray): n x m array
        centers (npt.NDArray): k x m array

    Returns:
        min_dist_to_centers (npt.NDArray): 1-D array of lenght n
    """
    n, m = data.shape
    k = centers.shape
    min_dist_to_centers = np.full(n, np.inf)
    for center in centers:
        dist = np.linalg.norm(data-center, axis=1)
        min_dist_to_centers = np.minimum(min_dist_to_centers, dist)
    return min_dist_to_centers


def farthestFirstTravel(data: npt.NDArray, k: int):
    """Given a set of data points, find k centers minimizing the maximum distance between these data points and centers.

    Args:
        data (NDArray): A numpy matrix, n x m, n points in m-dimensional space
        k (int): the number of centers
    """
    n, m = data.shape
    first = [0]  # [random.randrange(0, n)]
    centers_idx: set[int] = set(first)

    while len(centers_idx) < k:
        min_dist_to_centers = np.full(n, np.inf)

        for idx in centers_idx:
            center = data[idx, :]
            dist = np.linalg.norm(data-center, axis=1)
            min_dist_to_centers = np.minimum(min_dist_to_centers, dist)

        argmax = int(np.argmax(min_dist_to_centers))
        centers_idx.add(argmax)

    centers = data[np.array(list(centers_idx), dtype=np.int64), :]
    return centers


def kMeansInitializer(data: npt.NDArray, k: int):
    n, m = data.shape
    first = [random.randrange(0, n)]
    centers_idx: set[int] = set(first)

    while len(centers_idx) < k:
        idxs = np.array(list(centers_idx), dtype=np.int64)
        centers = data[idxs,:]
        min_dist_to_centers = minDistanceToCenters(data, centers)
        probs = min_dist_to_centers / sum(min_dist_to_centers)
        new_idx = np.random.choice(range(n), size=1, p=probs)
        centers_idx.add(int(new_idx))

    centers = data[np.array(list(centers_idx), dtype=np.int64), :]
    return centers


def LloydAlgorithm(data: npt.NDArray, k: int, 
                   init: Literal['kmeans_init', 'random', 'first_k'] = "kmeans_init"):
    """Implements Lloyd algorithm for k-means clustering

    Args:
        data (npt.NDArray): 
        k (int): number of clusters
        init (Literal['kmeans_init', 'random', 'first_k'], optional): How to initialize the algorithm. Defaults to "kmeans_init":
            - 'kmeans_init' uses the k-means initializer
            - 'random' randomly picks k points from the data to use as initial centers
            - 'first_k' picks the first k points in the data as initial centers

    Returns:
        centers (NDArray): an array with shape k x m 
    """
    n, m = data.shape

    # randomly allocate centers
    if init == 'kmeans_init':
        centers = kMeansInitializer(data, k)
    elif init == 'random':
        centers_idx = np.random.choice(n, size=k, replace=False)
        centers = data[centers_idx, :]
    else:
        centers = data[0:k,:]

    dist = np.inf
    while dist > 0:
        # assign each data point to a center (Centers to Clusters)
        assignment = argminDistanceToCenters(data, centers)

        # Move centers to the centers of clusters (Clusters to Centers)
        new_centers = np.zeros(centers.shape, dtype=centers.dtype)
        for i in range(k):
            cluster = data[assignment == i, :]
            center = np.mean(cluster, axis=0)
            new_centers[i, :] = center
        dist = np.linalg.norm(centers-new_centers)
        centers = new_centers

    return centers


def LloydAlgorithmMonteCarlo(data: npt.NDArray, k: int, iter: int):
    centers = LloydAlgorithm(data, k)
    distortion = squaredErrorDistortion(data, centers)

    for _ in range(iter-1):
        new_centers = LloydAlgorithm(data, k, 'kmeans_init')
        distortion2 = squaredErrorDistortion(data, new_centers)
        if distortion2 < distortion:
            centers = new_centers
            distortion = distortion2

    return centers


def softKMeansClustering(data: npt.NDArray, k: int, beta: float, iter: int = 100):
    """Implements the soft k-means clustering algorithm

    Args:
        data (npt.NDArray): _description_
        k (int): number of clusters
        beta (float): stiffness parameter
        iter (int, optional): Number of iterations. Defaults to 100.

    Returns:
        centers (npt.NDArray): k x m array of cluster centers
    """
    n, m = data.shape
    # select centers
    centers_idx = np.random.choice(n, size=k, replace=False)
    centers = data[centers_idx, :]

    for i in range(iter):
        # centers to soft clusters (E-step)
        # compute responsability matrix, k x n
        distances = distanceToCenters(data, centers)
        hidden_matrix = np.exp(-beta*distances).transpose()
        hidden_matrix = hidden_matrix / np.sum(hidden_matrix, axis=0)
        a = 0

        # soft clusters to centers (M-step)
        centers = (hidden_matrix @ data) / \
            np.sum(hidden_matrix, axis=1, keepdims=True)

    return centers
