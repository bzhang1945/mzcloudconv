# Use the official Python base image
FROM python:3.8-slim

# Set the working directory
WORKDIR /app

# Copy requirements.txt and install the dependencies
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application files
COPY . .

# Expose the port your app runs on
EXPOSE 8080

ENV SECRET_KEY 'b8e4fa1f06f3535a48f911bf611378a70c541211274f90c4de57687a29b94271'

# Start the application
CMD ["gunicorn", "-b", ":8080", "app:app", "--timeout", "120"]