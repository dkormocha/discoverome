FROM python:3.11.9

#system dependencies 
RUN apt-get update && apt-get install -y \
    build-essential \
    libpq-dev \
    curl \
    && rm -rf /var/lib/apt/lists/*

#install uv 
RUN pip install --no-cache-dir uv

#set workdir 
WORKDIR /app

#copy dependency files first (cache layer) 
COPY pyproject.toml uv.lock* /app/

#install dependencies using uv 
RUN uv sync --frozen

#copy application code 
COPY . .

#environment 
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    FLASK_APP=app/app.py

#expose port 
EXPOSE 6969

#run production server 
CMD ["uv", "run", "gunicorn", "-c", "app/gunicorn.conf.py", "app.app:app"]