"""
Author: Patryk Jarnot (2026)
"""
import sqlite3


class Cache:
    def __init__(self):
        self.db_path = "cache.sqlite"
        self.connection = None
        self.cursor = None
        self.init_database()

    def _tables_exist(self):
        pass

    def init_database(self):
        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()
        if not self._tables_exist():
            sql = """CREATE TABLE IF NOT EXISTS remoteitems (
                       id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                       category text NOT NULL,
                       param text NOT NULL,
                       value text NOT NULL);"""
            self.cursor.execute(sql)
            self.connection.commit()

    def store(self, category, param, value):
        param = param.replace(",", "_")
        self.cursor.execute('select value from remoteitems where category=? and param=?', (category, param))
        row = self.cursor.fetchone()
        if row is None:
            self.cursor.execute('insert into remoteitems (category, param, value) values (?, ?, ?)', (category, param, value))
        else:
            self.cursor.execute('update remoteitems set value=? where category=? and param=?', (value, category, param))
        self.connection.commit()

    def get(self, category, param):
        param = param.replace(",", "_")
        self.cursor.execute('select value from remoteitems where category=? and param=?', (category, param))
        row = self.cursor.fetchone()
        if row is None:
            return None
        else:
            return row[0]


