"""
Feature hierarchy for representing parent-child relationships in genomic features.
"""

from collections import defaultdict
from typing import List

from .record import GffRecord


class FeatureHierarchy:
    """
    Represents parent-child relationships between features.
    
    This class provides methods to navigate hierarchical relationships
    between features, particularly useful for GFF3 files with complex
    feature hierarchies.
    """
    
    def __init__(self):
        """Initialize an empty feature hierarchy."""
        self.parent_to_children = defaultdict(list)  # parent_id -> [child_records]
        self.child_to_parents = defaultdict(list)    # child_id -> [parent_records]
        self.id_to_record = {}                       # id -> record
    
    def add_record(self, record: GffRecord) -> None:
        """
        Add a record to the hierarchy.
        
        Args:
            record: The Record object to add
        """
        # Get record ID
        record_id = record.get_attribute('ID')
        if record_id:
            self.id_to_record[record_id] = record
        
        # Process parent-child relationships
        parents = record.get_attribute('Parent')
        if parents:
            # Convert to list if it's a single value
            if not isinstance(parents, list):
                parents = [parents]
                
            for parent_id in parents:
                self.child_to_parents[record_id].append(parent_id)
                self.parent_to_children[parent_id].append(record)
    
    def build_from_records(self, records: List[GffRecord]) -> None:
        """
        Build the hierarchy from a list of records.
        
        Args:
            records: List of Record objects
        """
        for record in records:
            self.add_record(record)
    
    def get_children(self, parent_id: str) -> List[GffRecord]:
        """
        Get all direct children of a feature.
        
        Args:
            parent_id: ID of the parent feature
            
        Returns:
            List of child Record objects
        """
        return self.parent_to_children.get(parent_id, [])
    
    def get_parents(self, child_id: str) -> List[GffRecord]:
        """
        Get all direct parents of a feature.
        
        Args:
            child_id: ID of the child feature
            
        Returns:
            List of parent Record objects
        """
        parent_ids = self.child_to_parents.get(child_id, [])
        return [self.id_to_record.get(pid) for pid in parent_ids if pid in self.id_to_record]
    
    def get_descendants(self, parent_id: str, max_depth: int | None = None) -> List[GffRecord]:
        """
        Get all descendants of a feature (children, grandchildren, etc.).
        
        Args:
            parent_id: ID of the parent feature
            max_depth: Maximum depth to traverse (None for unlimited)
            
        Returns:
            List of descendant Record objects
        """
        descendants = []
        visited = set()
        added_record_ids = set()
        
        def _traverse(pid, depth=0):
            if pid in visited or (max_depth is not None and depth >= max_depth):
                return
                
            visited.add(pid)
            children = self.get_children(pid)
            for child in children:
                child_id = child.get_attribute('ID')
                if child_id and child_id not in added_record_ids:
                    descendants.append(child)
                    added_record_ids.add(child_id)
            
            for child in children:
                child_id = child.get_attribute('ID')
                if child_id:
                    _traverse(child_id, depth + 1)
        
        _traverse(parent_id)
        return descendants
    
    def get_ancestors(self, child_id: str, max_depth: int | None = None) -> List[GffRecord]:
        """
        Get all ancestors of a feature (parents, grandparents, etc.).
        
        Args:
            child_id: ID of the child feature
            max_depth: Maximum depth to traverse (None for unlimited)
            
        Returns:
            List of ancestor Record objects
        """
        ancestors = []
        visited = set()
        
        def _traverse(cid, depth=0):
            if cid in visited or (max_depth is not None and depth >= max_depth):
                return
                
            visited.add(cid)
            parents = self.get_parents(cid)
            ancestors.extend(parents)
            
            for parent in parents:
                parent_id = parent.get_attribute('ID')
                if parent_id:
                    _traverse(parent_id, depth + 1)
        
        _traverse(child_id)
        return ancestors