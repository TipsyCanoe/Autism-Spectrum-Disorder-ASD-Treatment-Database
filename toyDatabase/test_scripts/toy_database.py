import faiss
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Generate sample data: 10 vectors with 3 dimensions
data = np.random.random((10, 3)).astype('float32')  # 10 vectors, 3 dimensions

# Build the FAISS index using L2 (Euclidean) distance
index = faiss.IndexFlatL2(3)  # 3D vectors
index.add(data)               # Add data to the index

# Create a random 3D query vector
query = np.random.random((1, 3)).astype('float32')
print("Query Vector:")
print(query)
print("-" * 30)

# Search for the top 3 nearest neighbors
distances, indices = index.search(query, k=3)  # Top 3 similar vectors
print("Nearest Matches (Indices):")
print(indices)
print("-" * 30)

# Display the most similar vectors with their distances
print("Most Similar Vectors:")
for idx, dist in zip(indices[0], distances[0]):
    print(f"Index: {idx}, Distance: {dist:.4f}, Vector: {data[idx]}")

# 3D Visualization
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot all data points
ax.scatter(data[:, 0], data[:, 1], data[:, 2], c='blue', marker='o', label='Data Points')

# Highlight the query vector
ax.scatter(query[0, 0], query[0, 1], query[0, 2], c='red', s=100, marker='^', label='Query Vector')

# Highlight the nearest neighbors
for idx in indices[0]:
    ax.scatter(data[idx, 0], data[idx, 1], data[idx, 2], c='green', s=80, marker='*', label='Nearest Neighbor')

# Draw lines from the query vector to its nearest neighbors
for idx in indices[0]:
    ax.plot([query[0, 0], data[idx, 0]],
            [query[0, 1], data[idx, 1]],
            [query[0, 2], data[idx, 2]],
            color='gray', linestyle='--')

# Labels and legend
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
ax.set_title('3D Visualization of Vector Similarity')

# To avoid duplicate labels in the legend
handles, labels = plt.gca().get_legend_handles_labels()
unique_labels = dict(zip(labels, handles))
plt.legend(unique_labels.values(), unique_labels.keys())

plt.show()
