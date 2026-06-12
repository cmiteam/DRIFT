#!/bin/bash
# Script to initialize a user with all base models

if [ -z "$1" ]; then
    echo "Usage: ./init_user.sh <username>"
    echo "Example: ./init_user.sh Rob"
    exit 1
fi

USERNAME=$1
USER_DIR="users/$USERNAME"

echo "Initializing user: $USERNAME"

# Create user directory structure
mkdir -p "$USER_DIR/models"

# Copy each base model
for MODEL in Default QuickTest Flood Creation OoA; do
    echo "  Creating model: $MODEL"

    MODEL_DIR="$USER_DIR/models/$MODEL"
    mkdir -p "$MODEL_DIR/config" "$MODEL_DIR/results"

    # Copy parameters from base model
    cp "static/basemodels/$MODEL/parameters.csv" "$MODEL_DIR/config/parameters.csv"

    # Get description from registry
    DESC=$(grep "^$MODEL," static/basemodels/registry.csv | cut -d',' -f3)

    # Create metadata
    cat > "$MODEL_DIR/metadata.json" <<EOF
{
  "name": "$MODEL",
  "description": "$DESC",
  "base_model": "$MODEL",
  "created": "$(date -Iseconds)",
  "last_modified": "$(date -Iseconds)",
  "author": "$USERNAME"
}
EOF
done

echo ""
echo "✓ User $USERNAME initialized successfully!"
echo "✓ Created 5 models: Default, QuickTest, Flood, Creation, OoA"
echo ""
echo "To run a simulation:"
echo "  ./drift.exe -username=$USERNAME -model=Default"
